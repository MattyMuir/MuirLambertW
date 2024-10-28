#include "w0.h"

#include <cstdint>
#include <limits>

#define EQUAL 0x0
#define LESS 0x11
#define GREATER 0x1E

#define SHUFFLE_INT(a, b, i) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), i))

// ========== Log Functions ==========
static inline __m256d Epi64ToPd(__m256i x)
{
    __m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    x = _mm256_permutevar8x32_epi32(x, perm);
    return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

static inline __m256d LogFast(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    // =================

    // Extract exponent
    __m256i punn = _mm256_castpd_si256(x);
    __m256d exp = Epi64ToPd(_mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023)));

    // Extract mantissa
    punn = _mm256_and_si256(punn, _mm256_set1_epi64x(0x3FFFFFFFFFFFFFFF));
    punn = _mm256_or_si256(punn, _mm256_set1_epi64x(0x3FF0000000000000));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Compute approximation
    static constexpr double P[] = {
        -1.4976283869142268,
        2.1229729992413895,
        -0.7362283025393341,
        0.11127972835311338
    };

    __m256d approx = _mm256_set1_pd(P[3]);
    for (size_t i = 0; i < 3; i++)
        approx = _mm256_fmadd_pd(approx, mantissa, _mm256_set1_pd(P[2 - i]));

    return _mm256_fmadd_pd(exp, ln2, approx);
}

static inline __m256d LogAccurate(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    __m256d one = _mm256_set1_pd(1.0);
    // =================

    // Extract exponent
    __m256d xd = _mm256_mul_pd(x, _mm256_set1_pd(1.0 / 0.75));
    __m256i dpunn = _mm256_castpd_si256(xd);
    __m256i exp = _mm256_sub_epi64(_mm256_srli_epi64(dpunn, 52), _mm256_set1_epi64x(1023));

    // Extract mantissa
    __m256i punn = _mm256_castpd_si256(x);
    punn = _mm256_sub_epi64(punn, _mm256_slli_epi64(exp, 52));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    __m256d t = _mm256_div_pd(_mm256_sub_pd(mantissa, one), _mm256_add_pd(mantissa, one));
    __m256d arg = _mm256_mul_pd(t, t);
    __m256d t3 = _mm256_mul_pd(t, arg);

    // Compute approximation
    static constexpr double P[] = {
        0.6666666666667778740063,
        0.399999999950799600689777,
        0.285714294746548025383248,
        0.222221366518767365905163,
        0.181863266251982985677316,
        0.152519917006351951593857,
        0.153487338491425068243146
    };

    __m256d approx = _mm256_set1_pd(P[6]);
    for (size_t i = 0; i < 6; i++)
        approx = _mm256_fmadd_pd(approx, arg, _mm256_set1_pd(P[5 - i]));
    approx = _mm256_fmadd_pd(t3, approx, _mm256_add_pd(t, t));

    return _mm256_fmadd_pd(Epi64ToPd(exp), ln2, approx);
}

// ========== General ==========
static inline __m256d Abs(__m256d x)
{
    __m256d signMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fff'ffff'ffff'ffff));
    return _mm256_and_pd(x, signMask);
}

static inline __m256d FirstApprox(__m256d x)
{
#if 1
    // === Constants ===
    __m256d e2 = _mm256_set1_pd(5.4365636569180904707);			// e * 2
    __m256d two = _mm256_set1_pd(2.0);							// 2
    // =================

    static constexpr double P[] = {
        -0.9998418255015216,
        0.9963519639934693,
        -0.31689383427598367,
        0.11927961836393368,
        -0.0389272342299362,
        0.009587149228046744,
        -0.0016619725895250465,
        0.000193387127759057,
        -1.426313511537818e-05,
        5.997840674383423e-07,
        -1.0923174871658654e-08
    };

    __m256d reta = _mm256_sqrt_pd(_mm256_fmadd_pd(x, e2, two));

    __m256d approx = _mm256_set1_pd(P[10]);
    for (size_t i = 0; i < 10; i++)
        approx = _mm256_fmadd_pd(approx, reta, _mm256_set1_pd(P[9 - i]));

    // Use approx = x for arguments near zero
    __m256d isNearZero = _mm256_cmp_pd(Abs(x), _mm256_set1_pd(1e-4), LESS);
    approx = _mm256_blendv_pd(approx, x, isNearZero);

    return approx;
#else
    static constexpr double P[] = {
        0,
        810.8417464502003,
        2372.4945068042407,
        1581.878759697207,
        224.10169908298582,
        4.420196050919542
    };

    static constexpr double Q[] = {
        810.8511630902697,
        3182.1361682278834,
        3545.937817132696,
        1172.6816431962184,
        92.41288592456523,
        1
    };

    __m256d numer = _mm256_set1_pd(P[5]);
    for (size_t i = 0; i < 5; i++)
        numer = _mm256_fmadd_pd(numer, x, _mm256_set1_pd(P[4 - i]));

    __m256d denom = _mm256_set1_pd(Q[5]);
    for (size_t i = 0; i < 5; i++)
        denom = _mm256_fmadd_pd(denom, x, _mm256_set1_pd(Q[4 - i]));

    __m256d approx = _mm256_div_pd(numer, denom);

    // Use approx = x for arguments near zero
    __m256d isNearZero = _mm256_cmp_pd(Abs(x), _mm256_set1_pd(1e-4), LESS);
    approx = _mm256_blendv_pd(approx, x, isNearZero);

    return approx;
#endif
}

static inline __m256d SecondApprox(__m256d x)
{
    static constexpr double P[] = {
        64393.137450661044568,
        43204.949550002405886,
        20295.724800471609342,
        453.37964270930132216,
        1.0
    };
    static constexpr double Q[] = {
        104344.40703457256313,
        22558.64516691800236,
        461.09954435682880103,
        0.9999372708768251572
    };

    __m256d logX = LogFast(x);

    __m256d numer = _mm256_set1_pd(P[4]);
    for (size_t i = 0; i < 4; i++)
        numer = _mm256_fmadd_pd(numer, logX, _mm256_set1_pd(P[3 - i]));

    __m256d denom = _mm256_set1_pd(Q[3]);
    for (size_t i = 0; i < 3; i++)
        denom = _mm256_fmadd_pd(denom, logX, _mm256_set1_pd(Q[2 - i]));

    __m256d approx = _mm256_div_pd(numer, denom);
    return approx;
}

static inline __m256d GeneralW0(__m256d x)
{
    __m256d isOver20 = _mm256_cmp_pd(x, _mm256_set1_pd(20.0), GREATER);
    uint32_t isOver20Mask = _mm256_movemask_pd(isOver20);

    __m256d w;
    switch (isOver20Mask)
    {
    case 0b0000:
        w = FirstApprox(x);
        break;
    case 0b1111:
        w = SecondApprox(x);
        break;
    default:
        w = _mm256_blendv_pd(FirstApprox(x), SecondApprox(x), isOver20);
    }

    // Constants
    __m256d c23 = _mm256_set1_pd(2.0 / 3.0);
    __m256d one = _mm256_set1_pd(1.0);

    // === Fritsch Iteration ===
    __m256d xov = _mm256_div_pd(x, w);
    __m256d zn = _mm256_sub_pd(LogAccurate(xov), w);

    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, c23, temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    return w;
}

// ========== Near Branch ==========
static inline __m256d AddEm(__m256d x)
{
    __m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
    __m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);

    return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

static inline __m256d NearBranchSeries(__m256d p)
{
    static constexpr double P[] = {
        -1,
        1.0,
        -0.3333333333333297,
        0.1527777777773187,
        -0.07962962960594742,
        0.04450231416408525,
        -0.025984703964189882,
        0.01563551165669095,
        -0.00961596052642643,
        0.006009378381399107,
        -0.0037902230498085155,
        0.002376673069259283,
        -0.0014296490725648598,
        0.000770653337014657,
        -0.00033482607204654753,
        9.991363954176918e-05,
        -1.4853323115219712e-05
    };

    // Evaluate polynomial using Horner's Method
    __m256d value = _mm256_set1_pd(P[16]);
    for (size_t i = 0; i < 16; i++)
        value = _mm256_fmadd_pd(value, p, _mm256_set1_pd(P[15 - i]));

    return value;
}

static inline __m256d NearBranchW0(__m256d x)
{
    static constexpr double s2e = 2.331643981597124;
    __m256d p = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(s2e));

    return NearBranchSeries(p);
}

// ========== Main Function ==========
__m256d MuirW0(__m256d x)
{
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.29), LESS);
    uint32_t nearBranchMask = _mm256_movemask_pd(isNearBranch);

    __m256d result;
    switch (nearBranchMask)
    {
    case 0b0000: [[likely]]
        result = GeneralW0(x);
        break;
    case 0b1111: [[unlikely]]
        result = NearBranchW0(x);
        break;
    default: [[unlikely]]
        result = _mm256_blendv_pd(GeneralW0(x), NearBranchW0(x), isNearBranch);
        break;
    }

    // Fix infinity
    __m256d infinity = _mm256_set1_pd(std::numeric_limits<double>::infinity());
    result = _mm256_blendv_pd(result, infinity, _mm256_cmp_pd(x, infinity, EQUAL));

    // Fix zero
    __m256d isZero = _mm256_cmp_pd(x, _mm256_setzero_pd(), EQUAL);
    result = _mm256_blendv_pd(result, x, isZero);

    return result;
}
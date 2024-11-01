#include "w0.h"

#include <cstdint>
#include <limits>

#define EQUAL 0x0
#define LESS 0x11
#define GREATER 0x1E

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
    static constexpr double P[] = {
        0,
        30.580056454638136,
        83.95836185597197,
        46.16620637664877,
        3.4636816277252214
    };

    static constexpr double Q[] = {
        30.578403642151667,
        114.49011569793561,
        114.80618615998705,
        28.635096582884064,
        1
    };

    __m256d numer = _mm256_set1_pd(P[4]);
    for (size_t i = 0; i < 4; i++)
        numer = _mm256_fmadd_pd(numer, x, _mm256_set1_pd(P[3 - i]));

    __m256d denom = _mm256_set1_pd(Q[4]);
    for (size_t i = 0; i < 4; i++)
        denom = _mm256_fmadd_pd(denom, x, _mm256_set1_pd(Q[3 - i]));

    __m256d approx = _mm256_div_pd(numer, denom);

    // Use approx = x for arguments near zero
    __m256d isNearZero = _mm256_cmp_pd(Abs(x), _mm256_set1_pd(1e-4), LESS);
    approx = _mm256_blendv_pd(approx, x, isNearZero);

    return approx;
}

static inline __m256d SecondApprox(__m256d x)
{
    static constexpr double P[] = {
        64312.7454007891,
        43264.12227598657,
        20243.65384336377,
        453.17656235798086,
        1.0000432316050645
    };
    static constexpr double Q[] = {
        104342.57917932322,
        22499.368605590193,
        460.93750724715477,
        1
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
    __m256d useLarge = _mm256_cmp_pd(x, _mm256_set1_pd(7.34), GREATER);
    uint32_t useLargeMask = _mm256_movemask_pd(useLarge);

    __m256d w;
    switch (useLargeMask)
    {
    case 0b0000:
        w = FirstApprox(x);
        break;
    case 0b1111:
        w = SecondApprox(x);
        break;
    default:
        w = _mm256_blendv_pd(FirstApprox(x), SecondApprox(x), useLarge);
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
        -1.00000000000000000000,
        0.99999999999998689937,
        -0.33333333333171155655,
        0.15277777769847986078,
        -0.07962962759798784818,
        0.04450228328389740917,
        -0.02598439214142129680,
        0.01563333375832150554,
        -0.00960508856297833703,
        0.00596982547465134492,
        -0.00368441824865070513,
        0.00216878673408957843,
        -0.00113330227139719539,
        0.00047252681627728467,
        -0.00013420111092875102,
        0.00001887878365359131,
    };

    // Evaluate polynomial using Horner's Method
    __m256d value = _mm256_set1_pd(P[15]);
    for (size_t i = 0; i < 15; i++)
        value = _mm256_fmadd_pd(value, p, _mm256_set1_pd(P[14 - i]));

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
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.28), LESS);
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
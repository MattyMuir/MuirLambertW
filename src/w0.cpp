#include <cstdint>
#include <limits>

#include <immintrin.h>

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
    __m256i exp = _mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023));

    // Extract mantissa
    punn = _mm256_sub_epi64(punn, _mm256_slli_epi64(exp, 52));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Polynomial approximation coefficients for algorithm 2, index 1, order 4
    static constexpr double P[] = {
        -1.74883843293730611990,
        2.84123605882110430443,
        -1.49089059863615602808,
        0.45671004396616393661,
        -0.05816743246573453929
    };

    __m256d approx = _mm256_set1_pd(P[4]);
    for (size_t i = 0; i < 4; i++)
        approx = _mm256_fmadd_pd(approx, mantissa, _mm256_set1_pd(P[3 - i]));

    return _mm256_fmadd_pd(Epi64ToPd(exp), ln2, approx);
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
static inline __m256d FirstApprox(__m256d x)
{
    // Rational approximation coefficients for algorithm 2, index 2, order 5/3
    static constexpr double P[] = {
        0,
        0.8115268950222906,
        2.4937647825894977,
        1.6258020048436266,
        0.08074281421946017,
        -0.0025394283830904394
    };
    static constexpr double Q[] = {
        0.8115268950222906,
        3.3053448229779,
        3.7148636123100385,
        1
    };

    __m256d numer = _mm256_set1_pd(P[5]);
    for (size_t i = 0; i < 5; i++)
        numer = _mm256_fmadd_pd(numer, x, _mm256_set1_pd(P[4 - i]));

    __m256d denom = _mm256_set1_pd(Q[3]);
    for (size_t i = 0; i < 3; i++)
        denom = _mm256_fmadd_pd(denom, x, _mm256_set1_pd(Q[2 - i]));

    __m256d approx = _mm256_div_pd(numer, denom);
    return approx;
}

static inline __m256d SecondApprox(__m256d x)
{
    // Rational approximation coefficients for algorithm 2, index 3, order 3/2
    static constexpr double P[] = {
        266.74662101711755,
        180.72015154289477,
        81.03709502548347,
        0.9987349881680496
    };
    static constexpr double Q[] = {
        438.5489337661638,
        87.11265453382124,
        1
    };

    __m256d logX = LogFast(x);

    __m256d numer = _mm256_set1_pd(P[3]);
    for (size_t i = 0; i < 3; i++)
        numer = _mm256_fmadd_pd(numer, logX, _mm256_set1_pd(P[2 - i]));

    __m256d denom = _mm256_set1_pd(Q[2]);
    for (size_t i = 0; i < 2; i++)
        denom = _mm256_fmadd_pd(denom, logX, _mm256_set1_pd(Q[1 - i]));

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
    // Polynomial approximation coefficients for algorithm 2, index 4, order 15
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

    return result;
}
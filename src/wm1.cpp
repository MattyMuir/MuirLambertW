#include <cstdint>
#include <limits>

#include <immintrin.h>

#define LESS 0x11
#define GREATER 0x1E

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
    __m256d denormScale = _mm256_set1_pd(4503599627370496.0);               // 2^52
    // =================

    // Fix for denormalized values
    x = _mm256_mul_pd(x, denormScale);

    // Extract exponent
    __m256i punn = _mm256_castpd_si256(x);
    __m256i exp = _mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023));

    // Extract mantissa
    punn = _mm256_sub_epi64(punn, _mm256_slli_epi64(exp, 52));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Polynomial approximation coefficients for algorithm 6, index 1, order 4
    static constexpr double P[] = {
        -37.772582545309206, //-1.7289291561920494e+00 - denormOffset (ln(2^52)),
        2.78901155791566960e+00,
        -1.44093748876198707e+00,
        4.36015488686681152e-01,
        -5.50266844824230939e-02
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

static inline __m256d Approx(__m256d x)
{
    // === Constants ===
    __m256d negOne = _mm256_set1_pd(-1.0);
    // =================

    // Rational approximation coefficients for algorithm 6, index 2, order 3/1
    static constexpr double P[] = {
        -3.836813614374928,
        -6.420188447784658,
        -3.950370152775716,
        -0.9985569509992126
    };
    static constexpr double Q[] = {
        3.8333830077075923,
        1.0
    };

    __m256d logX = LogFast(_mm256_sub_pd(_mm256_setzero_pd(), x));
    __m256d t = _mm256_sqrt_pd(_mm256_fmadd_pd(logX, negOne, negOne));

    __m256d numer = _mm256_set1_pd(P[3]);
    for (size_t i = 0; i < 3; i++)
        numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[2 - i]));

    __m256d denom = _mm256_add_pd(t, _mm256_set1_pd(Q[0]));

    return _mm256_div_pd(numer, denom);
}

static inline __m256d GeneralWm1(__m256d x)
{
    __m256d w = Approx(x);

    // === Constants ===
    __m256d c23 = _mm256_set1_pd(2.0 / 3.0);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d smallThreshold = _mm256_set1_pd(-1e-300);
    __m256d smallScale = _mm256_set1_pd(4611686018427387904.0);     // 2^62
    __m256d smallOffset = _mm256_set1_pd(42.975125194716609184);    // ln(2^62)
    // =================

    // === Fritsch Iteration ===
    // Compute zn
    __m256d isSmall = _mm256_cmp_pd(x, smallThreshold, GREATER);
    __m256d xScale = _mm256_blendv_pd(x, _mm256_mul_pd(x, smallScale), isSmall);
    __m256d xow = _mm256_div_pd(xScale, w);
    __m256d zn = LogAccurate(xow);
    zn = _mm256_blendv_pd(zn, _mm256_sub_pd(zn, smallOffset), isSmall);
    zn = _mm256_sub_pd(zn, w);

    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, c23, temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    return w;
}

static inline __m256d AddEm(__m256d x)
{
    __m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
    __m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);

    return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

static inline __m256d NearBranchSeries(__m256d p)
{
    // Polynomial approximation coefficients for algorithm 6, index 3, order 19
    static constexpr double P[] = {
        -1,
        -1.0000000000000016,
        -0.3333333333324043,
        -0.15277777790499153,
        -0.07962962215983704,
        -0.04450255318722647,
        -0.025980024318324864,
        -0.015697192504546754,
        -0.009049360327095492,
        -0.009821126316859557,
        0.01520694192144526,
        -0.0742821588262001,
        0.20519629823523558,
        -0.45489235269983486,
        0.7541893499759095,
        -0.9368341377306678,
        0.8408359779458122,
        -0.5187306631778058,
        0.19718629534561682,
        -0.03542966753367424
    };

    // Evaluate polynomial using Horner's Method
    __m256d value = _mm256_set1_pd(P[19]);
    for (size_t i = 0; i < 19; i++)
        value = _mm256_fmadd_pd(value, p, _mm256_set1_pd(P[18 - i]));

    return value;
}

static inline __m256d NearBranchWm1(__m256d x)
{
    static constexpr double s2e = 2.331643981597124;
    __m256d p = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(s2e));

    return NearBranchSeries(p);
}

__m256d MuirWm1(__m256d x)
{
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.27), LESS);
    uint32_t nearBranchMask = _mm256_movemask_pd(isNearBranch);

    __m256d result;
    switch (nearBranchMask)
    {
    case 0b0000: [[likely]]
        result = GeneralWm1(x);
        break;
    case 0b1111: [[unlikely]]
        result = NearBranchWm1(x);
        break;
    default: [[unlikely]]
        result = _mm256_blendv_pd(GeneralWm1(x), NearBranchWm1(x), isNearBranch);
        break;
    }

    return result;
}
#include "wm1.h"

#include <cstdint>
#include <limits>

#define LESS 0x11
#define GREATER 0x1E

static inline __m256d Epi64ToPd(__m256i x)
{
    __m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    x = _mm256_permutevar8x32_epi32(x, perm);
    return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

static __m256d LogFast(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    __m256d dblMin = _mm256_set1_pd(std::numeric_limits<double>::min());
    __m256d denormScale = _mm256_set1_pd(4503599627370496.0);               // 2^52
    __m256d denormOffset = _mm256_set1_pd(36.043653389117156090);           // ln(2^52)
    // =================

    // Fix for denormalized values
    __m256d isDenorm = _mm256_cmp_pd(x, dblMin, LESS);
    x = _mm256_blendv_pd(x, _mm256_mul_pd(x, denormScale), isDenorm);

    // Extract exponent
    __m256i punn = _mm256_castpd_si256(x);
    __m256d exp = Epi64ToPd(_mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023)));

    // Extract mantissa
    punn = _mm256_and_si256(punn, _mm256_set1_epi64x(0x3FFFFFFFFFFFFFFF));
    punn = _mm256_or_si256(punn, _mm256_set1_epi64x(0x3FF0000000000000));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Compute approximation
    static constexpr double P[] = {
        -1.7289291561920494e+00,
        2.78901155791566960e+00,
        -1.44093748876198707e+00,
        4.36015488686681152e-01,
        -5.50266844824230939e-02
    };

    __m256d approx = _mm256_set1_pd(P[4]);
    for (size_t i = 0; i < 4; i++)
        approx = _mm256_fmadd_pd(approx, mantissa, _mm256_set1_pd(P[3 - i]));

    // Fix for denormalized values
    approx = _mm256_blendv_pd(approx, _mm256_sub_pd(approx, denormOffset), isDenorm);

    return _mm256_fmadd_pd(exp, ln2, approx);
}

static __m256d LogAccurate(__m256d x)
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

__m256d Approx(__m256d x)
{
#if 0
    // === Constants ===
    __m256d zero = _mm256_setzero_pd();
    __m256d half = _mm256_set1_pd(0.5);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d negOne = _mm256_set1_pd(-1.0);
    __m256d six = _mm256_set1_pd(6.0);
    __m256d seven = _mm256_set1_pd(7.0);
    // =================

	// Compute t
    __m256d t = _mm256_sub_pd(negOne, LogAccurate(_mm256_sub_pd(zero, x)));

	// Compute x1
	static constexpr double P1[] = {
		3.00735538504012242e+00,
		-1.11079558616957819e-02,
		1.64300407704260418e-05,
		-2.87179186820882765e-08,
		2.93030279922565400e-11,
		-1.23057197181254966e-14
	};

    //3.00034140628619717*10^{+00},-5.15654510610554290*10^{-03},-1.13873366980070844*10^{-03},9.88231607008033285*10^{-05},-3.10110048923035139*10^{-06}
	static constexpr double P2[] = {
        3.00154719163439900e+00,
        -8.18302417096426201e-03,
        -1.59272314739640301e-04
	};

	__m256d x11 = _mm256_set1_pd(P1[5]);
	for (size_t i = 0; i < 5; i++)
		x11 = _mm256_fmadd_pd(x11, t, _mm256_set1_pd(P1[4 - i]));

    __m256d x12 = _mm256_set1_pd(P2[2]);
    for (size_t i = 0; i < 2; i++)
        x12 = _mm256_fmadd_pd(x12, t, _mm256_set1_pd(P2[1 - i]));

    __m256d x1 = _mm256_blendv_pd(x11, x12, _mm256_cmp_pd(t, _mm256_set1_pd(13.0), LESS));
    x1 = _mm256_div_pd(one, x1);

    __m256d approx = _mm256_sqrt_pd(_mm256_mul_pd(t, half));
    approx = _mm256_fmadd_pd(approx, x1, one);
    approx = _mm256_div_pd(six, approx);
    approx = _mm256_sub_pd(approx, t);
    approx = _mm256_sub_pd(approx, seven);

    return approx;
#else
    // Determine a
    __m256d xgt0 = _mm256_cmp_pd(x, _mm256_set1_pd(-0.29),      GREATER);
    __m256d xgt1 = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-11),     GREATER);
    __m256d xgt2 = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-96),     GREATER);

    __m256d a = _mm256_set1_pd(101.815);                    // [-1/e,       -0.29]
    a = _mm256_blendv_pd(a, _mm256_set1_pd(127), xgt0);     // [-0.29,      -1e-11]
    a = _mm256_blendv_pd(a, _mm256_set1_pd(181), xgt1);     // [-1e-11,     -1e-96]
    a = _mm256_blendv_pd(a, _mm256_set1_pd(317), xgt2);     // [-1e-96,     0]

    __m256d logX = LogFast(_mm256_sub_pd(_mm256_setzero_pd(), x));
    __m256d t = _mm256_sub_pd(_mm256_set1_pd(-1.0), logX);
    __m256d sqrtT = _mm256_sqrt_pd(t);
    __m256d approx = _mm256_fmadd_pd(a, sqrtT, _mm256_set1_pd(270));
    approx = _mm256_div_pd(t, approx);
    approx = _mm256_sub_pd(_mm256_set1_pd(1.0 / 3.0), approx);
    approx = _mm256_fmadd_pd(approx, sqrtT, _mm256_set1_pd(1.4142135623730950488));
    approx = _mm256_div_pd(sqrtT, approx);
    approx = _mm256_add_pd(approx, approx);
    approx = _mm256_sub_pd(logX, approx);

    return approx;
#endif
}

__m256d GeneralWm1(__m256d x)
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

__m256d AddEm(__m256d x)
{
    __m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
    __m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);

    return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

static __m256d NearBranchSeries(__m256d p)
{
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

static __m256d NearBranchWm1(__m256d x)
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
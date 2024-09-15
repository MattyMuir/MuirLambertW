#include "MuirpairWm1.h"

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
        -2.001981174757893,
        3.7491905722285463,
        -2.77383516599771,
        1.3494246360270206,
        -0.36414336579589207,
        0.04134449829706258
    };

    __m256d approx = _mm256_set1_pd(P[5]);
    for (size_t i = 0; i < 5; i++)
        approx = _mm256_fmadd_pd(approx, mantissa, _mm256_set1_pd(P[4 - i]));

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
#if 1
    // === Constants ===
    __m256d zero = _mm256_setzero_pd();
    __m256d half = _mm256_set1_pd(0.5);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d negOne = _mm256_set1_pd(-1.0);
    __m256d six = _mm256_set1_pd(6.0);
    __m256d seven = _mm256_set1_pd(7.0);
    // =================

	// Compute t
    __m256d t = _mm256_sub_pd(negOne, LogFast(_mm256_sub_pd(zero, x)));

	// Compute x1
	static constexpr double P1[] = {
		3.00735538504012242e+00,
		-1.11079558616957819e-02,
		1.64300407704260418e-05,
		-2.87179186820882765e-08,
		2.93030279922565400e-11,
		-1.23057197181254966e-14
	};

	static constexpr double P2[] = {
		3.00142773717038159e+00,
		-8.00996540237069126e-03,
		-1.80448428405709059e-04
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
    __m256d xQuiteSmall = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-11), GREATER);
    __m256d xVerySmall = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-96), GREATER);
    __m256d a = _mm256_blendv_pd(_mm256_set1_pd(127), _mm256_set1_pd(181), xQuiteSmall);
    a = _mm256_blendv_pd(a, _mm256_set1_pd(317), xVerySmall);

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

__m256d MuirpairWm1(__m256d x)
{
    __m256d w = Approx(x);
    //return w;

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
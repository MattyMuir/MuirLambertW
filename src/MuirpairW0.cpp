#include "MuirpairW0.h"

#include "mysleef.h"

#define GREATER 0x1E

#define SHUFFLE_INT(a, b, i) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), i))

static __m256d LogApprox(__m256d x)
{
    __m256i punn = _mm256_castpd_si256(x);
    __m128i punnTop = _mm256_extractf128_si256(punn, 1);
    __m128i compressed = _mm256_castsi256_si128(SHUFFLE_INT(punn, _mm256_castsi128_si256(punnTop), 221));
    compressed = _mm_sub_epi32(compressed, _mm_set1_epi32(1072632447));
    return _mm256_mul_pd(_mm256_cvtepi32_pd(compressed), _mm256_set1_pd(6.610368362777016e-7));
}

__m256d FirstApprox(__m256d x)
{
	// === Constants ===
	__m256d e2 = _mm256_set1_pd(5.4365636569180904707);			// e * 2
	__m256d one = _mm256_set1_pd(1.0);
	__m256d two = _mm256_set1_pd(2.0);							// 2
	__m256d three = _mm256_set1_pd(3.0);						// 3
	__m256d b = _mm256_set1_pd(1.09556884765625);				// Barry et al. p.164
	__m256d a = _mm256_set1_pd(4.612634277343749);				// Barry et al. p.164
	__m256d rt2m1 = _mm256_set1_pd(0.41421356237309504880);		// sqrt(2) - 1
	__m256d m = _mm256_set1_pd(0.29289321881345247560);			// (2*sqrt(2)-3)/(sqrt(2)-2)
	// =================

	__m256d reta = _mm256_sqrt_pd(_mm256_fmadd_pd(x, e2, two));
	__m256d n2 = _mm256_sqrt_pd(_mm256_sqrt_pd(_mm256_add_pd(reta, b)));
	n2 = _mm256_mul_pd(n2, a);
	__m256d n1 = _mm256_fmadd_pd(n2, m, rt2m1);
	__m256d d = _mm256_div_pd(_mm256_mul_pd(n1, reta), _mm256_add_pd(n2, reta));
	__m256d approx = _mm256_add_pd(d, three);
	approx = _mm256_div_pd(reta, approx);
	approx = _mm256_add_pd(approx, one);
	approx = _mm256_div_pd(reta, approx);
	approx = _mm256_sub_pd(approx, one);

	return approx;
}

__m256d SecondApprox(__m256d x)
{
    // === Constants ===
    __m256d y = _mm256_set1_pd(-5.70115661621093750e+00);
    __m256d offset = _mm256_set1_pd(-0.132);
    // =================

    static constexpr double P[] = {
       6.42275660145116698e+00,
       1.33047964073367945e+00,
       6.72008923401652816e-02,
       1.16444069958125895e-03,
       7.06966760237470501e-06,
       5.48974896149039165e-09,
       -7.00379652018853621e-11,
       -1.89247635913659556e-13,
       -1.55898770790170598e-16,
       -4.06109208815303157e-20,
       -2.21552699006496737e-24,
    };
    static constexpr double Q[] = {
       1.00000000000000000e+00,
       3.34498588416632854e-01,
       2.51519862456384983e-02,
       6.81223810622416254e-04,
       7.94450897106903537e-06,
       4.30675039872881342e-08,
       1.10667669458467617e-10,
       1.31012240694192289e-13,
       6.53282047177727125e-17,
       1.11775518708172009e-20,
       3.78250395617836059e-25,
    };

    __m256d logX = LogApprox(_mm256_add_pd(x, offset));

    __m256d numer = _mm256_set1_pd(P[8]);
    __m256d denom = _mm256_set1_pd(Q[8]);
    for (size_t i = 0; i < 8; i++)
    {
        numer = _mm256_fmadd_pd(numer, logX, _mm256_set1_pd(P[7 - i]));
        denom = _mm256_fmadd_pd(denom, logX, _mm256_set1_pd(Q[7 - i]));
    }

    __m256d approx = _mm256_div_pd(numer, denom);
    approx = _mm256_add_pd(_mm256_add_pd(approx, y), logX);
    return approx;
}

__m256d MuirpairW0(__m256d x)
{
    __m256d approx1 = FirstApprox(x);
    __m256d approx2 = SecondApprox(x);
    __m256d isOver20 = _mm256_cmp_pd(x, _mm256_set1_pd(20.0), GREATER);
    __m256d w = _mm256_blendv_pd(approx1, approx2, isOver20);;

    // Constants
    __m256d c23 = _mm256_set1_pd(2.0 / 3.0);
    __m256d one = _mm256_set1_pd(1.0);

    // === Fritsch Iteration ===
    __m256d xov = _mm256_div_pd(x, w);
    __m256d zn = _mm256_sub_pd(Sleef_logd4_u35avx2(xov), w);

    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, c23, temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    return w;
}
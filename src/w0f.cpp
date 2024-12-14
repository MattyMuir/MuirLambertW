#include <cfloat>
#include <cstdint>

#include <immintrin.h>

#define LESS 0x11

struct double8
{
    __m256d l, h;
};

static inline double8 Set1(double x)
{
	return { _mm256_set1_pd(x), _mm256_set1_pd(x) };
}

static inline double8 fmadd(double8 x, double8 y, double z)
{
    __m256d zp = _mm256_set1_pd(z);
    return { _mm256_fmadd_pd(x.l, y.l, zp), _mm256_fmadd_pd(x.h, y.h, zp) };
}

static inline double8 div(double8 x, double8 y)
{
	return { _mm256_div_pd(x.l, y.l), _mm256_div_pd(x.h, y.h) };
}

static inline double8 FloatToDouble(__m256 p)
{
    __m256d l = _mm256_cvtps_pd(_mm256_castps256_ps128(p));
    __m256d h = _mm256_cvtps_pd(_mm256_extractf128_ps(p, 1));
    return { l, h };
}

static inline __m256 DoubleToFloat(double8 p)
{
    return _mm256_set_m128(_mm256_cvtpd_ps(p.h), _mm256_cvtpd_ps(p.l));
}

static inline __m256 LogAccurate(__m256 x)
{
	// === Constants ===
	__m256 one = _mm256_set1_ps(1.0f);
	__m256 ln2 = _mm256_set1_ps(0.693147180559945286226764f);
	// =================

	// Extract exponent
	__m256 xd = _mm256_mul_ps(x, _mm256_set1_ps(1.0f / 0.75f));
	__m256i dpunn = _mm256_castps_si256(xd);
	__m256i exp = _mm256_sub_epi32(_mm256_srli_epi32(dpunn, 23), _mm256_set1_epi32(127));

	// Extract mantissa
	__m256i punn = _mm256_castps_si256(x);
	punn = _mm256_sub_epi32(punn, _mm256_slli_epi32(exp, 23));
	__m256 mantissa = _mm256_castsi256_ps(punn);

	// Evaluate series
	__m256 t = _mm256_div_ps(_mm256_sub_ps(mantissa, one), _mm256_add_ps(mantissa, one));
	__m256 t2 = _mm256_mul_ps(t, t);

	static constexpr float P[] = {
		2.0f,
		0.666666686534881591796875f,
		0.400005877017974853515625f,
		0.28518211841583251953125f,
		0.2392828464508056640625f
	};

	__m256 approx = _mm256_set1_ps(P[4]);
	for (size_t i = 0; i < 4; i++)
		approx = _mm256_fmadd_ps(approx, t2, _mm256_set1_ps(P[3 - i]));

	return _mm256_fmadd_ps(t, approx, _mm256_mul_ps(ln2, _mm256_cvtepi32_ps(exp)));
}

static inline __m256 AddEm(__m256 x)
{
	const __m256 emHigh = _mm256_set1_ps(0.36787945f);
	const __m256 emLow = _mm256_set1_ps(-9.149756e-09f);
	return _mm256_add_ps(_mm256_add_ps(x, emHigh), emLow);
}

static inline __m256 NearBranchW0(__m256 x)
{
	__m256 rt2e = _mm256_set1_ps(2.331644f);

	static constexpr float P[] = {
		-0.9999999781289544,
		0.9999966080647236,
		-0.33324531164727067,
		0.15189891604646868,
		-0.07530393941472714,
		0.03290035332102544,
		-0.008369773627101843
	};

	__m256 p = _mm256_mul_ps(_mm256_sqrt_ps(AddEm(x)), rt2e);

	__m256 res = _mm256_set1_ps(P[6]);
	for (size_t i = 0; i < 6; i++)
		res = _mm256_fmadd_ps(res, p, _mm256_set1_ps(P[5 - i]));

	return res;
}

static inline __m256 FirstApprox(__m256 x)
{
	static constexpr double P[] = {
		0,
		152.91909285006932,
		1024.8955779609623,
		2455.092024342446,
		2528.414753137949,
		1066.824695440063,
		148.54000057530493,
		4.0539143590175835
	};

	static constexpr double Q[] = {
		152.91909349763372,
		1177.8147226529784,
		3403.5283401899756,
		4573.002542975176,
		2878.913736848507,
		761.147965396105,
		65.5632478027466,
		1
	};

	double8 xd = FloatToDouble(x);

	double8 numer = Set1(P[7]);
	for (size_t i = 0; i < 7; i++)
		numer = fmadd(numer, xd, P[6 - i]);

	double8 denom = Set1(Q[7]);
	for (size_t i = 0; i < 7; i++)
		denom = fmadd(denom, xd, Q[6 - i]);

	return DoubleToFloat(div(numer, denom));
}

static inline __m256 SecondApprox(__m256 x)
{
	static constexpr double P[] = {
		245182.20097823755,
		280243.5212428723,
		142843.813324628,
		40353.72076097795,
		5776.914448840662,
		184.83613670644033,
		0.9984483567344636
	};

	static constexpr double Q[] = {
		432788.26007218857,
		216948.13159273885,
		58081.26591912717,
		6594.751582203545,
		191.21022696372594,
		1
	};

	double8 t = FloatToDouble(LogAccurate(x));

	double8 numer = Set1(P[6]);
	for (size_t i = 0; i < 6; i++)
		numer = fmadd(numer, t, P[5 - i]);

	double8 denom = Set1(Q[5]);
	for (size_t i = 0; i < 5; i++)
		denom = fmadd(denom, t, Q[4 - i]);

	return DoubleToFloat(div(numer, denom));
}

static inline __m256 GeneralW0(__m256 x)
{
	__m256 useFirst = _mm256_cmp_ps(x, _mm256_set1_ps(6.9035267829895019531f), LESS);
	uint32_t useFirstMask = _mm256_movemask_ps(useFirst);

	__m256 result;
	switch (useFirstMask)
	{
	case 0b00000000: [[likely]]
		result = SecondApprox(x);
		break;
	case 0b11111111: [[unlikely]]
		result = FirstApprox(x);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(SecondApprox(x), FirstApprox(x), useFirst);
		break;
	}

	return result;
}

__m256 MuirW0(__m256 x)
{
	__m256 isNearBranch = _mm256_cmp_ps(x, _mm256_set1_ps(-0.3f), LESS);
	uint32_t nearBranchMask = _mm256_movemask_ps(isNearBranch);

	__m256 result;
	switch (nearBranchMask)
	{
	case 0b00000000: [[likely]]
		result = GeneralW0(x);
		break;
	case 0b11111111: [[unlikely]]
		result = NearBranchW0(x);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(GeneralW0(x), NearBranchW0(x), isNearBranch);
		break;
	}

	return result;
}
#include <cfloat>
#include <cmath>
#include <cstdint>

#include <immintrin.h>

#define EQUAL 0x0
#define LESS 0x11

static inline __m256 LogApprox(__m256 x)
{
	// === Constants ===
	__m256i off = _mm256_set1_epi32(-1064866805);
	__m256 scale = _mm256_set1_ps(8.262958405176314e-8f);
	// =================

	__m256i punn = _mm256_castps_si256(x);
	punn = _mm256_add_epi32(punn, off);
	return _mm256_mul_ps(_mm256_cvtepi32_ps(punn), scale);
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

static inline __m256 FirstApprox(__m256 x)
{
	static constexpr float P[] = {
		0.0f,
		1.0173324233f,
		1.62516706451f,
	};

	static constexpr float Q[] = {
		1.0173324233f,
		2.61423411949f,
		1.0f
	};

	__m256 numer = _mm256_set1_ps(P[2]);
	for (size_t i = 0; i < 2; i++)
		numer = _mm256_fmadd_ps(numer, x, _mm256_set1_ps(P[1 - i]));

	__m256 denom = _mm256_set1_ps(Q[2]);
	for (size_t i = 0; i < 2; i++)
		denom = _mm256_fmadd_ps(denom, x, _mm256_set1_ps(Q[1 - i]));

	return _mm256_div_ps(numer, denom);
}

static inline __m256 SecondApprox(__m256 x)
{
	static constexpr float P[] = {
		3.69555171276f,
		2.77242714985f,
		0.987426086605f
	};

	static constexpr float Q = 6.45416961763;

	__m256 t = LogApprox(x);

	__m256 numer = _mm256_set1_ps(P[2]);
	for (size_t i = 0; i < 2; i++)
		numer = _mm256_fmadd_ps(numer, t, _mm256_set1_ps(P[1 - i]));

	__m256 denom = _mm256_add_ps(t, _mm256_set1_ps(Q));

	return _mm256_div_ps(numer, denom);
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
		-0.999999966739467f,
		0.9999951977820332f,
		-0.3332171041026461f,
		0.15169332645744893f,
		-0.07462906014376834f,
		0.03188270977312735f,
		-0.0077961216354129935f
	};

	__m256 p = _mm256_mul_ps(_mm256_sqrt_ps(AddEm(x)), rt2e);

	__m256 res = _mm256_set1_ps(P[6]);
	for (size_t i = 0; i < 6; i++)
		res = _mm256_fmadd_ps(res, p, _mm256_set1_ps(P[5 - i]));

	return res;
}

static inline __m256 GeneralW0(__m256 x)
{
	// === Constants ===
	__m256 c23 = _mm256_set1_ps(2.0 / 3.0);
	__m256 one = _mm256_set1_ps(1.0);
	// =================

	__m256 useFirst = _mm256_cmp_ps(x, _mm256_set1_ps(1.6487212707), LESS);
	__m256 w;
	switch (_mm256_movemask_ps(useFirst))
	{
	case 0b00000000:
		w = SecondApprox(x);
		break;
	case 0b11111111:
		w = FirstApprox(x);
		break;
	default:
		w = _mm256_blendv_ps(SecondApprox(x), FirstApprox(x), useFirst);
	}

	// === Fritsch Iteration ===
	__m256 xov = _mm256_div_ps(x, w);
	__m256 zn = _mm256_sub_ps(LogAccurate(xov), w);

	__m256 temp = _mm256_add_ps(w, one);
	__m256 temp2 = _mm256_fmadd_ps(zn, c23, temp);
	temp2 = _mm256_mul_ps(temp, temp2);
	temp2 = _mm256_add_ps(temp2, temp2);
	__m256 temp3 = _mm256_sub_ps(temp2, _mm256_add_ps(zn, zn));
	__m256 temp4 = _mm256_mul_ps(_mm256_div_ps(zn, temp), _mm256_sub_ps(temp2, zn));
	w = _mm256_mul_ps(w, _mm256_add_ps(_mm256_div_ps(temp4, temp3), one));

	return w;
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

	// Fix infinity
	__m256 infinity = _mm256_set1_ps(INFINITY);
	result = _mm256_blendv_ps(result, infinity, _mm256_cmp_ps(x, infinity, EQUAL));

	return result;
}
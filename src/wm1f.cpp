#include <cstdint>
#include <cfloat>

#include <immintrin.h>

#define LESS 0x11
#define BLEND_INT(a, b, mask) _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), mask))

static inline __m256 LogAccurate(__m256 x)
{
	// === Constants ===
	__m256 one = _mm256_set1_ps(1.0f);
	__m256 denormScale = _mm256_set1_ps(1.8446744E19f);
	__m256i denormOffset = _mm256_set1_epi32(64);
	__m256 ln2 = _mm256_set1_ps(0.693147180559945286226764f);
	// =================

	// Scale x for denormalized values
	__m256 isDenorm = _mm256_cmp_ps(x, _mm256_set1_ps(FLT_MIN), LESS);
	x = _mm256_blendv_ps(x, _mm256_mul_ps(x, denormScale), isDenorm);

	// Extract exponent
	__m256 xd = _mm256_mul_ps(x, _mm256_set1_ps(1.0f / 0.75f));
	__m256i dpunn = _mm256_castps_si256(xd);
	__m256i exp = _mm256_sub_epi32(_mm256_srli_epi32(dpunn, 23), _mm256_set1_epi32(127));

	// Extract mantissa
	__m256i punn = _mm256_castps_si256(x);
	punn = _mm256_sub_epi32(punn, _mm256_slli_epi32(exp, 23));
	__m256 mantissa = _mm256_castsi256_ps(punn);

	// Fix exponent for denormalized values
	exp = BLEND_INT(exp, _mm256_sub_epi32(exp, denormOffset), isDenorm);

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

static inline __m256 NearBranchWm1(__m256 x)
{
	// === Constants ===
	__m256 rt2e = _mm256_set1_ps(2.331644f);
	// =================

	// Compute p
	__m256 p = _mm256_mul_ps(_mm256_sqrt_ps(AddEm(x)), rt2e);

	// Polynomial approximation coefficients for algorithm 8, index 3, order 8
	static constexpr float P[] = {
		-1.0000000212255025,
		-0.9999950257094332,
		-0.3335221441109622,
		-0.1500624689844515,
		-0.09883211897476347,
		0.029915601909860912,
		-0.1889854031081707,
		0.177910291677026,
		-0.11266484471740043
	};

	__m256 res = _mm256_set1_ps(P[8]);
	for (size_t i = 0; i < 8; i++)
		res = _mm256_fmadd_ps(res, p, _mm256_set1_ps(P[7 - i]));

	return res;
}

static inline __m256 FirstApprox(__m256 t)
{
	// Polynomial approximation coefficients for algorithm 8, index 1, order 7
	static constexpr float P[] = {
		-0.9999947770248371,
		-1.4142802274833077,
		-0.6663468462971003,
		-0.07936575583288796,
		0.016000547102702947,
		-0.0024352660263118973,
		0.00024032743757630923,
		-1.1250272634544634e-05
	};

	__m256 res = _mm256_set1_ps(P[7]);
	for (size_t i = 0; i < 7; i++)
		res = _mm256_fmadd_ps(res, t, _mm256_set1_ps(P[6 - i]));

	return res;
}

static inline __m256 SecondApprox(__m256 t)
{
	// Polynomial approximation coefficients for algorithm 8, index 2, order 7
	static constexpr float P[] = {
		-1.0918049,
		-1.2658587,
		-0.77344185,
		-0.033979293,
		0.0036363874,
		-0.0002548953,
		1.0409415e-05,
		-1.8701131e-07
	};

	__m256 res = _mm256_set1_ps(P[7]);
	for (size_t i = 0; i < 7; i++)
		res = _mm256_fmadd_ps(res, t, _mm256_set1_ps(P[6 - i]));

	return res;
}

static inline __m256 GeneralWm1(__m256 x)
{
	// === Constants ===
	__m256 negOne = _mm256_set1_ps(-1.0f);
	// =================

	// Compute t
	__m256 t = LogAccurate(_mm256_sub_ps(_mm256_setzero_ps(), x));
	t = _mm256_sqrt_ps(_mm256_sub_ps(negOne, t));

	__m256 useFirst = _mm256_cmp_ps(x, _mm256_set1_ps(-0.00000224905596703), LESS);
	uint32_t useFirstMask = _mm256_movemask_ps(useFirst);

	__m256 result;
	switch (useFirstMask)
	{
	case 0b00000000:
		result = SecondApprox(t);
		break;
	case 0b11111111:
		result = FirstApprox(t);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(SecondApprox(t), FirstApprox(t), useFirst);
		break;
	}

	return result;
}

__m256 MuirWm1(__m256 x)
{
	__m256 isNearBranch = _mm256_cmp_ps(x, _mm256_set1_ps(-0.307278738601), LESS);
	uint32_t nearBranchMask = _mm256_movemask_ps(isNearBranch);

	__m256 result;
	switch (nearBranchMask)
	{
	case 0b00000000: [[likely]]
		result = GeneralWm1(x);
		break;
	case 0b11111111: [[unlikely]]
		result = NearBranchWm1(x);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(GeneralWm1(x), NearBranchWm1(x), isNearBranch);
		break;
	}

	return result;
}
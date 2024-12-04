#include "wm1f.h"

#include <cstdint>
#include <cfloat>

#define LESS 0x11
#define BLEND_INT(a, b, mask) _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), mask))

__m256 LogAccurate(__m256 x)
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

__m256 AddEm(__m256 x)
{
	const __m256 emHigh = _mm256_set1_ps(0.36787945f);
	const __m256 emLow = _mm256_set1_ps(-9.149756e-09f);
	return _mm256_add_ps(_mm256_add_ps(x, emHigh), emLow);
}

__m256 NearBranchWm1(__m256 x)
{
	// === Constants ===
	__m256 rt2e = _mm256_set1_ps(2.331644f);
	// =================

	// Compute p
	__m256 p = _mm256_mul_ps(_mm256_sqrt_ps(AddEm(x)), rt2e);

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

__m256 GeneralWm1(__m256 x)
{
	// === Constants ===
	__m256 negTwo = _mm256_set1_ps(-2.0f);
	// =================

	// Compute t
	__m256 t = LogAccurate(_mm256_sub_ps(_mm256_setzero_ps(), x));
	t = _mm256_sqrt_ps(_mm256_fmadd_ps(t, negTwo, negTwo));

	static constexpr float P[] = {
		-0.9999492423798767,
		-1.0002532686721437,
		-0.33279221980879625,
		-0.028436955137531802,
		0.004219360823001216,
		-0.0005080457276522693,
		4.65793262714409e-05,
		-2.9923469210392722e-06,
		1.1273454917460403e-07,
		-6.441939908321432e-10,
		-1.562062365956019e-10,
		7.014581753864945e-12,
		-1.0236004967997851e-13
	};

	__m256 res = _mm256_set1_ps(P[12]);
	for (size_t i = 0; i < 12; i++)
		res = _mm256_fmadd_ps(res, t, _mm256_set1_ps(P[11 - i]));

	return res;
}

__m256 MuirWm1(__m256 x)
{
	__m256 isNearBranch = _mm256_cmp_ps(x, _mm256_set1_ps(-0.277689970954), LESS);
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
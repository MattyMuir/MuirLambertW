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

// [EM, -0.307278738601]
static inline __m256 NearBranchWm1(__m256 x)
{
	static constexpr float P[] = {
		-0.9999999820862636899409037,-2.331653425770667928726464,-1.811384580994698926256565,-1.96222368838886514852528,-1.961393601522052749815649,-6.23046445702523402424888,9.180235188603716180047407,-31.30984139409068550600034
	};

	__m256 t = _mm256_sqrt_ps(AddEm(x));

	__m256 res = _mm256_set1_ps(P[7]);
	for (size_t i = 0; i < 7; i++)
		res = _mm256_fmadd_ps(res, t, _mm256_set1_ps(P[6 - i]));

	return res;
}

// [-0.307278738601, -0.00000224905596703]
static inline __m256 FirstApprox(__m256 lx)
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

	// Compute t
	__m256 negOne = _mm256_set1_ps(-1.0f);
	__m256 t = _mm256_sqrt_ps(_mm256_sub_ps(negOne, lx));

	__m256 res = _mm256_set1_ps(P[7]);
	for (size_t i = 0; i < 7; i++)
		res = _mm256_fmadd_ps(res, t, _mm256_set1_ps(P[6 - i]));

	return res;
}

// [-0.00000224905596703, 0]
static inline __m256 SecondApprox(__m256 lx)
{
	static constexpr float P[] = {
		2103.5388959664333427,-8754.041759120255355,677.4414443307139027,-6.307429022866525761,0.0015828923769514411498
	};
	static constexpr float Q[] = {
		-9738.437009819998806785,3751.2359976224974621098,-165.07802916229234665566,1.
	};

	__m256 numer = _mm256_set1_ps(P[4]);
	for (size_t i = 0; i < 4; i++)
		numer = _mm256_fmadd_ps(numer, lx, _mm256_set1_ps(P[3 - i]));

	__m256 denom = _mm256_set1_ps(Q[3]);
	for (size_t i = 0; i < 3; i++)
		denom = _mm256_fmadd_ps(denom, lx, _mm256_set1_ps(Q[2 - i]));

	return _mm256_add_ps(_mm256_div_ps(numer, denom), lx);
}

static inline __m256 GeneralWm1(__m256 x)
{
	// Compute lx
	__m256 lx = LogAccurate(_mm256_sub_ps(_mm256_setzero_ps(), x));

	__m256 useFirst = _mm256_cmp_ps(x, _mm256_set1_ps(-0.00000224905596703), LESS);
	uint32_t useFirstMask = _mm256_movemask_ps(useFirst);

	__m256 result;
	switch (useFirstMask)
	{
	case 0b00000000:
		result = SecondApprox(lx);
		break;
	case 0b11111111:
		result = FirstApprox(lx);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(SecondApprox(lx), FirstApprox(lx), useFirst);
		break;
	}

	return result;
}

__m256 MuirWm1v2(__m256 x)
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
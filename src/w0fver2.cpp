#include <cfloat>
#include <cstdint>

#include <immintrin.h>

#define LESS 0x11

static inline __m256 LogApprox(__m256 x)
{
	// === Constants ===
	__m256 ln2 = _mm256_set1_ps(0.69314718055994530942f);
	// =================

	// Extract exponent
	__m256i punn = _mm256_castps_si256(x);
	__m256i exp = _mm256_sub_epi32(_mm256_srli_epi32(punn, 23), _mm256_set1_epi32(127));

	// Extract mantissa
	punn = _mm256_sub_epi32(punn, _mm256_slli_epi32(exp, 23));
	__m256 mantissa = _mm256_castsi256_ps(punn);

	// Compute approximation
	static constexpr float P[] = {
		-1.8143835131817594,
		3.031865632321167,
		-1.6945291478989106,
		0.5513886009540612,
		-0.07434157219377215
	};

	__m256 approx = _mm256_set1_ps(P[4]);
	for (size_t i = 0; i < 4; i++)
		approx = _mm256_fmadd_ps(approx, mantissa, _mm256_set1_ps(P[3 - i]));

	return _mm256_fmadd_ps(_mm256_cvtepi32_ps(exp), ln2, approx);
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

static inline __m256 Pow(__m256i q)
{
	return _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_add_epi32(q, _mm256_set1_epi32(0x7f)), 23));
}

static inline __m256 LoadExp(__m256 x, __m256i e)
{
	return _mm256_mul_ps(x, Pow(e));
}

static inline __m256 ExpAccurate(__m256 x)
{
	// === Constants ===
	__m256 rln2 = _mm256_set1_ps(1.44269504089f);
	__m256 l2u = _mm256_set1_ps(-0.693145751953125f);
	__m256 l2l = _mm256_set1_ps(-1.428606765330187045e-06f);
	__m256 one = _mm256_set1_ps(1.0f);
	// =================

	__m256i q = _mm256_cvtps_epi32(_mm256_mul_ps(x, rln2));
	__m256 s, u;

	s = _mm256_fmadd_ps(_mm256_cvtepi32_ps(q), l2u, x);
	s = _mm256_fmadd_ps(_mm256_cvtepi32_ps(q), l2l, s);

	u = _mm256_set1_ps(0.000198527617612853646278381);
	u = _mm256_fmadd_ps(u, s, _mm256_set1_ps(0.00139304355252534151077271));
	u = _mm256_fmadd_ps(u, s, _mm256_set1_ps(0.00833336077630519866943359));
	u = _mm256_fmadd_ps(u, s, _mm256_set1_ps(0.0416664853692054748535156));
	u = _mm256_fmadd_ps(u, s, _mm256_set1_ps(0.166666671633720397949219));
	u = _mm256_fmadd_ps(u, s, _mm256_set1_ps(0.5));

	u = _mm256_add_ps(one, _mm256_fmadd_ps(_mm256_mul_ps(s, s), u, s));

	u = LoadExp(u, q);

	return u;
}

static inline __m256 Abs(__m256 x)
{
	__m256i signMask = _mm256_set1_epi32(0x7fffffff);
	return _mm256_castsi256_ps(_mm256_and_si256(_mm256_castps_si256(x), signMask));
}

static inline __m256 Approx(__m256 x)
{
	__m256 linearThreshold = _mm256_set1_ps(0.01f);

	static constexpr float P[] = {
		0.0f,
		7.42689907805,
		13.2300488085,
		4.41802326989,
		0.977801471478,
		0.0000848767774065
	};

	static constexpr float Q[] = {
		7.42739098061,
		16.9250195294,
		7.91421158475,
		1.0
	};

	__m256 t = LogApprox(_mm256_add_ps(x, _mm256_set1_ps(1.0f)));

	__m256 numer = _mm256_set1_ps(P[5]);
	for (size_t i = 0; i < 5; i++)
		numer = _mm256_fmadd_ps(numer, t, _mm256_set1_ps(P[4 - i]));

	__m256 denom = _mm256_set1_ps(Q[3]);
	for (size_t i = 0; i < 3; i++)
		denom = _mm256_fmadd_ps(denom, t, _mm256_set1_ps(Q[2 - i]));

	__m256 approx = _mm256_div_ps(numer, denom);
	__m256 isLinear = _mm256_cmp_ps(Abs(x), linearThreshold, LESS);
	return _mm256_blendv_ps(approx, x, isLinear);
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
	__m256 one = _mm256_set1_ps(1.0f);
	__m256 two = _mm256_set1_ps(2.0f);
	// =================

	__m256 w = Approx(x);

	// === Halley Iteration ===
	__m256 ew = ExpAccurate(w);
	__m256 wew = _mm256_mul_ps(w, ew);
	__m256 wewx = _mm256_fmsub_ps(w, ew, x);
	__m256 w1 = _mm256_add_ps(w, one);
	__m256 w2 = _mm256_add_ps(w, two);
	__m256 res = _mm256_div_ps(_mm256_mul_ps(w2, wewx), _mm256_fmadd_ps(w, two, two));
	res = _mm256_fmsub_ps(ew, w1, res);
	res = _mm256_sub_ps(w, _mm256_div_ps(wewx, res));

	return res;
}

__m256 MuirW0v2(__m256 x)
{
	__m256 isNearBranch = _mm256_cmp_ps(x, _mm256_set1_ps(-0.29f), LESS);
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
#include <cmath>
#include <cstdint>
#include <immintrin.h>

#define EQUAL 0x0
#define LESS 0x11

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

// [EM, -0.205466397497]
static inline __m256 Approx1(__m256 t)
{
	static constexpr float P[] = {
		-7.5383643613783287804368,0.75093710445194423096623,16.0853573089437500092647,5.2230201577831681638114
	};
	static constexpr float Q[] = {
		7.5383643693344303338354,16.8258426675850213904358,9.4856801812012282150965,1.
	};

	__m256 numer = _mm256_set1_ps(P[3]);
	for (size_t i = 0; i < 3; i++)
		numer = _mm256_fmadd_ps(numer, t, _mm256_set1_ps(P[2 - i]));

	__m256 denom = _mm256_set1_ps(Q[3]);
	for (size_t i = 0; i < 3; i++)
		denom = _mm256_fmadd_ps(denom, t, _mm256_set1_ps(Q[2 - i]));

	return _mm256_div_ps(numer, denom);
}

// [-0.205466397497, 25.5337217474]
static inline __m256 Approx2(__m256 x, __m256 t)
{
	static constexpr float P[] = {
		-10973.510261722770450288,23229.03330548828134504,21732.6338059439559095,-1758.10647708696266683,-924.4396303454805579442,-2.
	};
	static constexpr float Q[] = {
		-29830.310685986906640648,-29465.753207256893110428,-5436.1389751146039623073,-73.338483763659317060332,1.
	};

	__m256 numer = _mm256_set1_ps(P[5]);
	for (size_t i = 0; i < 5; i++)
		numer = _mm256_fmadd_ps(numer, t, _mm256_set1_ps(P[4 - i]));

	__m256 denom = _mm256_set1_ps(Q[4]);
	for (size_t i = 0; i < 4; i++)
		denom = _mm256_fmadd_ps(denom, t, _mm256_set1_ps(Q[3 - i]));

	__m256 two = _mm256_set1_ps(2.0f);
	__m256 res = _mm256_div_ps(numer, denom);
	res = _mm256_fmadd_ps(t, two, res);
	return _mm256_div_ps(x, res);
}

// [25.5337217474, FLT_MAX]
static inline __m256 Approx3(__m256 x)
{
	static constexpr float P[] = {
		-97.60169098646942404,135.11605369354352048,-104.86329655432149538,71.1427070578416942,-4.035664801369620418,-0.13376645763002767975,0.00160027854040215705
	};
	static constexpr float Q[] = {
		-105.29972040791455897365,24.441774729384367267568,-30.459795510157873770576,-9.237885712512041305418,1.
	};

	__m256 lx = LogAccurate(x);
	__m256 t = _mm256_sqrt_ps(lx);

	__m256 numer = _mm256_set1_ps(P[6]);
	for (size_t i = 0; i < 6; i++)
		numer = _mm256_fmadd_ps(numer, t, _mm256_set1_ps(P[5 - i]));

	__m256 denom = _mm256_set1_ps(Q[4]);
	for (size_t i = 0; i < 4; i++)
		denom = _mm256_fmadd_ps(denom, t, _mm256_set1_ps(Q[3 - i]));

	return _mm256_add_ps(_mm256_div_ps(numer, denom), lx);
}

static inline __m256 NearOriginW0(__m256 x)
{
	__m256 t = _mm256_sqrt_ps(AddEm(x));

	__m256 isNearBranch = _mm256_cmp_ps(x, _mm256_set1_ps(-0.205466397497f), LESS);
	uint32_t nearBranchMask = _mm256_movemask_ps(isNearBranch);

	__m256 result;
	switch (nearBranchMask)
	{
	case 0b00000000: [[likely]]
		result = Approx2(x, t);
		break;
	case 0b11111111: [[unlikely]]
		result = Approx1(t);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(Approx2(x, t), Approx1(t), isNearBranch);
		break;
	}

	return result;
}

__m256 MuirW0(__m256 x)
{
	__m256 isNearOrigin = _mm256_cmp_ps(x, _mm256_set1_ps(25.5337217474f), LESS);
	uint32_t nearOriginMask = _mm256_movemask_ps(isNearOrigin);

	__m256 result;
	switch (nearOriginMask)
	{
	case 0b00000000:
		result = Approx3(x);
		break;
	case 0b11111111:
		result = NearOriginW0(x);
		break;
	default: [[unlikely]]
		result = _mm256_blendv_ps(Approx3(x), NearOriginW0(x), isNearOrigin);
		break;
	}

	// Fix infinity
	__m256 infinity = _mm256_set1_ps(INFINITY);
	result = _mm256_blendv_ps(result, infinity, _mm256_cmp_ps(x, infinity, EQUAL));

	return result;
}
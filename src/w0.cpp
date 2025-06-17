#include <cstdint>
#include <cmath>
#include <immintrin.h>

#define EQUAL 0x0
#define LESS 0x11
#define GREATER 0x1E

static inline __m256d Epi64ToPd(__m256i x)
{
	__m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
	x = _mm256_permutevar8x32_epi32(x, perm);
	return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

static inline __m256d LogAccurate(__m256d x)
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

static inline __m256d AddEm(__m256d x)
{
	__m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
	__m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);
	return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

// [EM, 7]
static inline __m256d NearOriginW0(__m256d x)
{
	static constexpr double P[] = {
		-0.00712055882855771169474,-0.6703572146960279620638,-2.24448365927303885572217,-2.51494773946536910201475,-0.85425231935584144352953,0.23719595996884174434909,0.191567663961857499505996,0.0314843726863854913180816,0.001300214864537770506208227,4.615781083048832162273736e-7
	};
	static constexpr double Q[] = {
		1.,3.949282680450614688385034,5.820473886693447585668792,4.000932412205590469209318,1.313258423915893827346632,0.1886092985291080619818581,0.009118162792373766308391816,0.00004198356440537148146361118,-3.077187388699221441515824e-7
	};

	__m256d t = _mm256_sqrt_pd(AddEm(x));

	__m256d numer = _mm256_set1_pd(P[9]);
	for (size_t i = 0; i < 9; i++)
		numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[8 - i]));

	__m256d denom = _mm256_set1_pd(Q[8]);
	for (size_t i = 0; i < 8; i++)
		denom = _mm256_fmadd_pd(denom, t, _mm256_set1_pd(Q[7 - i]));

	__m256d res = _mm256_div_pd(numer, denom);
	res = _mm256_add_pd(res, _mm256_set1_pd(0.375));
	res = _mm256_fmadd_pd(t, _mm256_set1_pd(1.5), res);
	return _mm256_div_pd(x, res);
}

// [7, 882046.451763]
static inline __m256d Approx2(__m256d lx)
{
	static constexpr double P[] = {
		1.989018967036558218269332,0.966994579900128877003661,0.289256463612459600820318,0.0217256029729141586290329,-0.0057684373966958188633333,-0.00203029254008971618759777,-0.000124832647437063624918513,-5.647695181765254876187e-8,3.629635022952366205546e-8,3.18108853738725602516e-11
	};
	static constexpr double Q[] = {
		1.,0.8069819561459024544431381,0.3672706764235715757396559,0.09953764313174786210693154,0.01678825816345157269586683,0.001479072475985500870448731,0.00003743216433752187452558595,-4.153601532168373105288119e-7,-8.975049417438828872303809e-9
	};

	__m256d numer = _mm256_set1_pd(P[9]);
	for (size_t i = 0; i < 9; i++)
		numer = _mm256_fmadd_pd(numer, lx, _mm256_set1_pd(P[8 - i]));

	__m256d denom = _mm256_set1_pd(Q[8]);
	for (size_t i = 0; i < 8; i++)
		denom = _mm256_fmadd_pd(denom, lx, _mm256_set1_pd(Q[7 - i]));

	__m256d res = _mm256_div_pd(numer, denom);
	return _mm256_add_pd(_mm256_sub_pd(res, _mm256_set1_pd(1.421875)), lx);
}

// [882046.451763, DBL_MAX]
static inline __m256d Approx3(__m256d lx)
{
	static constexpr double P[] = {
		-1.00533881015610352e6,2.3911456970362179e6,-30169.087551103194,-133375.547205836751,615877.19212552896,62470.99168047287,-6826.5341075058529,-558.42585300721034,-7.3294963867575946,-0.0065103995006925356,7.8588470804253449e-6
	};
	static constexpr double Q[] = {
		-58554.722745651911714408,101831.38386680150401492,345944.1074680087980566,-114053.5644056862896046,95619.16859950479625675,49842.145571153604896002,4985.9037608534460772842,143.80477919346177699662,1.
	};

	__m256d t = _mm256_sqrt_pd(lx);

	__m256d numer = _mm256_set1_pd(P[10]);
	for (size_t i = 0; i < 10; i++)
		numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[9 - i]));

	__m256d denom = _mm256_set1_pd(Q[8]);
	for (size_t i = 0; i < 8; i++)
		denom = _mm256_fmadd_pd(denom, t, _mm256_set1_pd(Q[7 - i]));

	__m256d res = _mm256_div_pd(numer, denom);
	return _mm256_add_pd(_mm256_sub_pd(res, _mm256_set1_pd(4.4892578125)), lx);
}

static inline __m256d GeneralW0(__m256d x)
{
	__m256d lx = LogAccurate(x);

	__m256d useLarge = _mm256_cmp_pd(x, _mm256_set1_pd(882046.451763), GREATER);
	uint32_t useLargeMask = _mm256_movemask_pd(useLarge);

	__m256d w;
	switch (useLargeMask)
	{
	case 0b0000:
		w = Approx2(lx);
		break;
	case 0b1111:
		w = Approx3(lx);
		break;
	default: [[unlikely]]
		w = _mm256_blendv_pd(Approx2(lx), Approx3(lx), useLarge);
	}

	return w;
}

__m256d MuirW0(__m256d x)
{
	__m256d isNearOrigin = _mm256_cmp_pd(x, _mm256_set1_pd(7), LESS);
	uint32_t nearOriginMask = _mm256_movemask_pd(isNearOrigin);

	__m256d w;
	switch (nearOriginMask)
	{
	case 0b0000:
		w = GeneralW0(x);
		break;
	case 0b1111:
		w = NearOriginW0(x);
		break;
	default: [[unlikely]]
		w = _mm256_blendv_pd(GeneralW0(x), NearOriginW0(x), isNearOrigin);
		break;
	}

	// Fix infinity
	__m256d infinity = _mm256_set1_pd(INFINITY);
	w = _mm256_blendv_pd(w, infinity, _mm256_cmp_pd(x, infinity, EQUAL));

	return w;
}
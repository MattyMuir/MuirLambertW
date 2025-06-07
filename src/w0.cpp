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

// [7, 162754.791419]
static inline __m256d Approx2(__m256d lx)
{
	static constexpr double P[] = {
		2.067143967078050320967014,1.069806179885573906650058,0.337809405204211931957211,0.0356729188661552557646515,-0.0038535621323684274778649,-0.00199257943881031994971427,-0.000158094881578911012337975,-2.7318652677779251602255e-6,-3.081471407150281460288e-9,4.953176349882239750154106e-12
	};
	static constexpr double Q[] = {
		1.,0.8262191934468344673637481,0.3828164400267764276286399,0.1066359739504771743003248,0.01872978639074896942897624,0.001812960554120384996829358,0.00006863753172739671577627016,6.73673778448057724709718e-7
	};

	__m256d numer = _mm256_set1_pd(P[9]);
	for (size_t i = 0; i < 9; i++)
		numer = _mm256_fmadd_pd(numer, lx, _mm256_set1_pd(P[8 - i]));

	__m256d denom = _mm256_set1_pd(Q[7]);
	for (size_t i = 0; i < 7; i++)
		denom = _mm256_fmadd_pd(denom, lx, _mm256_set1_pd(Q[6 - i]));

	__m256d res = _mm256_div_pd(numer, denom);
	return _mm256_add_pd(_mm256_sub_pd(res, _mm256_set1_pd(1.5)), lx);
}

// [162754.791419, DBL_MAX]
static inline __m256d Approx3(__m256d lx)
{
	static constexpr double P[] = {
		5.678562481658061702610738,-1.038809946424879152819072,5.034891416154293393930052,-0.265064170881860308933986,0.764486826236326682525158,0.515393261091868449856043,0.0271008407556604663865665,-0.0032768185098398277196103,-0.00019106602268372965491118,-2.017323947244474399957e-6,-1.52330111348782608258e-9,1.5738544705959299029e-12
	};
	static constexpr double Q[] = {
		1.,-0.05182544086903837325584622,0.711483346193023155259492,0.2813455070046693803259518,0.03654754123742773693916102,0.131185068186327629236442,0.03232285548397857649537308,0.002245093918805030411525548,0.00004976917202276323164806426,2.799409091991322856983097e-7
	};

	__m256d t = _mm256_sqrt_pd(lx);

	__m256d numer = _mm256_set1_pd(P[11]);
	for (size_t i = 0; i < 11; i++)
		numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[10 - i]));

	__m256d denom = _mm256_set1_pd(Q[9]);
	for (size_t i = 0; i < 9; i++)
		denom = _mm256_fmadd_pd(denom, t, _mm256_set1_pd(Q[8 - i]));

	__m256d res = _mm256_div_pd(numer, denom);
	return _mm256_add_pd(_mm256_sub_pd(res, _mm256_set1_pd(5)), lx);
}

static inline __m256d GeneralW0(__m256d x)
{
	__m256d lx = LogAccurate(x);

	__m256d useLarge = _mm256_cmp_pd(x, _mm256_set1_pd(162754.791419), GREATER);
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
#include "ReferenceLambertW.h"

#include <iostream>
#include <format>
#include <bit>

#include <mpfr.h>

int64_t GetExp(double x)
{
	uint64_t punn = std::bit_cast<uint64_t>(x);
	int64_t exp = (int64_t)((punn >> 52) & 0x7ff) - 1023;
	return exp;
}

Interval Bisection(double x, mpfr_t low, mpfr_t high, bool increasing, mpfr_prec_t precision)
{
	mpfr_t m;
	mpfr_init2(m, 53);

	mpfr_t yLow, yHigh;
	mpfr_init2(yLow, precision);
	mpfr_init2(yHigh, precision);
	for (;;)
	{
		// m = (low + high) / 2
		mpfr_add(m, low, high, MPFR_RNDN);
		mpfr_div_2ui(m, m, 1, MPFR_RNDN);

		if (mpfr_equal_p(m, low) || mpfr_equal_p(m, high))
			break; // Bracket cannot be narrowed any further

		// Compute yLow
		mpfr_rnd_t rnd = (mpfr_cmp_ui(m, 0) > 0) ? MPFR_RNDD : MPFR_RNDU;
		mpfr_exp(yLow, m, rnd);
		mpfr_mul(yLow, yLow, m, MPFR_RNDD);
		mpfr_sub_d(yLow, yLow, x, MPFR_RNDD);

		// Compute yHigh
		rnd = (mpfr_cmp_ui(m, 0) > 0) ? MPFR_RNDU : MPFR_RNDD;
		mpfr_exp(yHigh, m, rnd);
		mpfr_mul(yHigh, yHigh, m, MPFR_RNDU);
		mpfr_sub_d(yHigh, yHigh, x, MPFR_RNDU);

		int lowCmp = mpfr_cmp_ui(yLow, 0);
		int highCmp = mpfr_cmp_ui(yHigh, 0);

		if (lowCmp >= 0 && highCmp >= 0)
			mpfr_set(increasing ? high : low, m, MPFR_RNDN);
		if (lowCmp <= 0 && highCmp <= 0)
			mpfr_set(increasing ? low : high, m, MPFR_RNDN);

		if (lowCmp < 0 && highCmp > 0)
		{
			std::cerr << std::format("Error, ambiguous sign : {}\n", x);
			throw;
		}
	}

	double lowD = mpfr_get_d(low, MPFR_RNDD);
	double highD = mpfr_get_d(high, MPFR_RNDU);

	mpfr_clear(m);
	mpfr_clear(yLow);
	mpfr_clear(yHigh);

	return { lowD, highD };
}

Interval ReferenceLambertW0(double x)
{
	// Constants
	static constexpr double Infinity = std::numeric_limits<double>::infinity();
	static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

	// 1.5700924586837747e-16
	static constexpr mpfr_prec_t W0Precision = 107;

	// Edge cases
	if (x < EM_UP)
		return { NaN, NaN };
	if (x == Infinity)
		return { std::numeric_limits<double>::max(), Infinity };

	// Initialize variables
	mpfr_t xMpfr, low, high;
	mpfr_init2(xMpfr, 53);
	mpfr_init2(low, 53);
	mpfr_init2(high, 53);

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	// high = ln(x + 1)
	mpfr_log1p(high, xMpfr, MPFR_RNDU);

	if (x > 14)
	{
		// low = (E(x + 2) - E(ln2 * E(x + 2) + 1) - 1) * ln2
		static constexpr double ln2 = 0.69314718055994530942;
		int64_t e1 = GetExp(x + 2);
		int64_t e2 = GetExp(e1 * ln2 + 1.0);
		double bound = (e1 - e2 - 1) * ln2;
		mpfr_set_d(low, bound, MPFR_RNDN);
	}
	else if (x >= 0)
	{
		// low = x / (x + 1)
		mpfr_t xp1;
		mpfr_init2(xp1, 53);
		mpfr_add_ui(xp1, xMpfr, 1, MPFR_RNDU);
		mpfr_div(low, xMpfr, xp1, MPFR_RNDD);
		mpfr_clear(xp1);
	}
	else
	{
		// low = x * (1 - x * 8)
		mpfr_t x8;
		mpfr_init2(x8, 53);
		mpfr_mul_2ui(x8, xMpfr, 3, MPFR_RNDN);
		mpfr_ui_sub(low, 1, x8, MPFR_RNDU);
		mpfr_mul(low, low, xMpfr, MPFR_RNDD);
		mpfr_clear(x8);
		
		// Clamp low above -1
		if (mpfr_cmp_si(low, -1) < 0)
			mpfr_set_si(low, -1, MPFR_RNDN);
	}

	// === Bisection ===
	auto ret = Bisection(x, low, high, true, W0Precision);

	mpfr_clear(xMpfr);
	mpfr_clear(low);
	mpfr_clear(high);

	return ret;
}

Interval ReferenceLambertWm1(double x)
{
	// -0.3678794304802449
	static constexpr mpfr_prec_t Wm1Precision = 82;

	// Edge cases
	static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
	if (x < EM_UP || x >= 0)
		return { NaN, NaN };

	// Initialize variables
	mpfr_t xMpfr, low, high;
	mpfr_init2(xMpfr, 53);
	mpfr_init2(low, 53);
	mpfr_init2(high, 53);

	// Convert x to mpfr
	mpfr_set_d(xMpfr, x, MPFR_RNDN);

	// === Compute Bracket ===
	mpfr_t uDown;
	mpfr_init2(uDown, 53);
	mpfr_neg(uDown, xMpfr, MPFR_RNDN);
	mpfr_log(uDown, uDown, MPFR_RNDU);
	mpfr_neg(uDown, uDown, MPFR_RNDN);
	mpfr_sub_ui(uDown, uDown, 1, MPFR_RNDD);

	mpfr_t uUp;
	mpfr_init2(uUp, 53);
	mpfr_neg(uUp, xMpfr, MPFR_RNDN);
	mpfr_log(uUp, uUp, MPFR_RNDD);
	mpfr_neg(uUp, uUp, MPFR_RNDN);
	mpfr_sub_ui(uUp, uUp, 1, MPFR_RNDU);

	// low = -1 - (sqrt(u * 2) + u);
	mpfr_mul_2ui(low, uUp, 1, MPFR_RNDU);
	mpfr_sqrt(low, low, MPFR_RNDU);
	mpfr_add(low, low, uUp, MPFR_RNDU);
	mpfr_si_sub(low, -1, low, MPFR_RNDD);

	// high = -1 - (sqrt(u * 2) + u * 2 / 3)
	mpfr_mul_2ui(high, uDown, 1, MPFR_RNDD);
	mpfr_sqrt(high, high, MPFR_RNDD);
	mpfr_mul_2ui(uDown, uDown, 1, MPFR_RNDD);
	mpfr_div_ui(uDown, uDown, 3, MPFR_RNDD);
	mpfr_add(high, high, uDown, MPFR_RNDD);
	mpfr_si_sub(high, -1, high, MPFR_RNDU);

	mpfr_clear(uDown);
	mpfr_clear(uUp);

	// === Bisection ===
	auto ret = Bisection(x, low, high, false, Wm1Precision);

	mpfr_clear(xMpfr);
	mpfr_clear(low);
	mpfr_clear(high);

	return ret;
}
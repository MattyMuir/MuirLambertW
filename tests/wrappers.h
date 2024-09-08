#pragma once
#include <numeric>

#include "types.h"

static constexpr double Infinity = std::numeric_limits<double>::infinity();

template <MpfrFunction1D Func>
Interval MakeBounded(double x)
{
	// Initialize xMp
	mpfr_t xMp;
	mpfr_init2(xMp, 53);
	mpfr_set_d(xMp, x, MPFR_RNDN);

	// Evaluate function
	mpfr_t yMp;
	mpfr_init2(yMp, 53);
	int ternary = Func(yMp, xMp, MPFR_RNDD);
	double y = mpfr_get_d(yMp, MPFR_RNDD);

	// Free memory
	mpfr_clear(xMp);
	mpfr_clear(yMp);

	if (ternary == 0)
		return { y, y };
	if (ternary < 0)
		return { y, std::nextafter(y, INFINITY) };

	throw;
}

template <SimdFunction1D Func>
double MakeSerial(double x)
{
	__m256d p = _mm256_set1_pd(x);
	__m256d res = Func(p);
	return res[0];
}
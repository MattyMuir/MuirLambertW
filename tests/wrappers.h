#pragma once
#include <numeric>

#include "types.h"

static constexpr double Infinity = std::numeric_limits<double>::infinity();

template <SimdFunction1D Func>
double MakeSerial(double x)
{
	__m256d p = _mm256_set1_pd(x);
	__m256d res = Func(p);
	return res[0];
}
#pragma once
#include <immintrin.h>

#include <functional>

#include <ReferenceW.h>

// Function types
using Function1D = double(*)(double);
using BoundedFunction1D = Interval(*)(double);
using MpfrFunction1D = int(*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
using SimdFunction1D = __m256d(*)(__m256d);

using RandomFunction = std::function<double()>;
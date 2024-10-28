#pragma once
#include <immintrin.h>

#include <functional>

#include <ReferenceW.h>

// Function types
using Function1D = double(*)(double);
using BoundedFunction1D = Interval(*)(double);
using SimdFunction1D = __m256d(*)(__m256d);

using RandomFunction = std::function<double()>;
#pragma once
#include <immintrin.h>

#include <functional>

#include <mpfr.h>

#include <ReferenceLambertW.h>

// === Function types ===
// float
using Function1Df = float(*)(float);
using BoundedFunction1Df = Intervalf(*)(float);
using SimdFunction1Df = __m256(*)(__m256);
using RandomFunctionf = std::function<float()>;

// double
using Function1D = double(*)(double);
using BoundedFunction1D = Interval(*)(double);
using SimdFunction1D = __m256d(*)(__m256d);
using RandomFunction = std::function<double()>;

// MPFR
using MpfrFunction1D = int(*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
#include <iostream>
#include <vector>
#include <random>
#include <format>
#include <chrono>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"
#include "constants.h"

#include "MuirLambertW.h"
#include <ReferenceLambertW.h>

Intervalf ReferenceW0f(float x)
{
	static ReferenceWf evaluator;
	return evaluator.W0(x);
}

Intervalf ReferenceWm1f(float x)
{
	static ReferenceWf evaluator;
	return evaluator.Wm1(x);
}

Interval ReferenceW0(double x)
{
	static ReferenceW evaluator;
	return evaluator.W0(x);
}

Interval ReferenceWm1(double x)
{
	static ReferenceW evaluator;
	return evaluator.Wm1(x);
}

float ExpMapWm1(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP / (1 + exp(x));
}

float MuirWm1MadeSerial(float x)
{
	__m256 p = _mm256_set1_ps(x);
	int res = _mm_extract_ps(_mm256_castps256_ps128(MuirWm1(p)), 0);
	return std::bit_cast<float>(res);
}

int main()
{
#if 1
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx<float> dist{ EM_UPf, 0, false };

	MaxULPRounded(ReferenceWm1f, MuirWm1MadeSerial, [&]() { return dist(gen); });
#endif

	ULPHistogramSigned(ReferenceWm1f, MuirWm1MadeSerial, -15, 50, 0.2, ExpMapWm1);
}
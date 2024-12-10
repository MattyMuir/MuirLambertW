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

double ExpMapWm1(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	return EM_UP / (1 + exp(x));
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

double MuirW0MadeSerial(double x)
{
	__m256d p = _mm256_set1_pd(x);
	int64_t res = _mm256_extract_epi64(_mm256_castpd_si256(MuirW0(p)), 0);
	return std::bit_cast<double>(res);
}

int main()
{
#if 1
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx<float> dist{ EM_UPf, 0, false };

	MaxULPRounded(ReferenceWm1f, MuirWm1MadeSerial, [&]() { return dist(gen); });
#else
	float x = EM_UPf;

	for (size_t i = 0; x < 0; i++)
	{
		auto wApprox = MuirWm1MadeSerial(x);
		auto wExact = ReferenceWm1f(x);

		uint32_t err = ULPDistance(wApprox, wExact);
		if (err > 4)
		{
			std::cout << std::format("Error: {}\n", x);
			break;
		}

		x = std::nextafter(x, INFINITY);
		if (i % 100'000 == 0)
			std::cout << x << '\n';
	}
#endif
}
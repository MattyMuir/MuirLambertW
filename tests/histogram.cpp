#include "histogram.h"

#include <iostream>
#include <format>
#include <random>

#include "ulp.h"

void ULPHistogram(BoundedFunction1D boundedFunc, Function1D approxFunc, double min, double max, double step)
{
	static std::mt19937_64 gen{ std::random_device{}() };

	for (double low = min; low < max; low += step)
	{
		double high = low + step;
		std::uniform_real_distribution<double> dist{ low, high };

		uint64_t err = MaxULPRounded(boundedFunc, approxFunc, [&]() { return dist(gen); }, 10'000);
		std::cout << std::format("{:.5f} - {:.5f}: {}\n", low, high, err);
	}
}
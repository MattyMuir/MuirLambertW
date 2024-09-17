#include <vector>
#include <random>
#include <format>
#include <chrono>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"

int main()
{
	/*
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ EM_UP, -0.25 };

	double worst = -1;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = ReferenceLambertW0(x);

		double approx = MuirW0(x);

		uint64_t err = std::max(ULPDistance(approx, inf), ULPDistance(approx, sup));

		if (err <= 4)
			continue;

		if (x > worst)
		{
			worst = x;
			dist = std::uniform_real_distribution<double>{ worst, -0.25 };
			std::cout << std::format("{}\n", worst);
		}
	}
	*/

	static std::mt19937_64 gen{ std::random_device{}() };
	static std::uniform_real_distribution<double> dist{ EM_UP, 0 };

	size_t NumData = 10'000;

	std::vector<double> data;
	for (size_t i = 0; i < NumData; i++)
		data.push_back(dist(gen));

	auto start = std::chrono::steady_clock::now();
	double _ = 0.0;
	for (double d : data)
		_ += ReferenceLambertW0(d).inf;
	auto end = std::chrono::steady_clock::now();

	std::cout << _ << '\n';
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start) << '\n';
}
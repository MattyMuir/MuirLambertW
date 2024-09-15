#include <random>
#include <format>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"

double BoostLambertW0(double x)
{
	return boost::math::lambert_w0(x);
}

int main()
{
	/*
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ EM_UP, 0 };

	double worst = -1;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = ReferenceLambertWM1(x);

		double approx = BarryLambertWM1(x);

		uint64_t dist = std::max(ULPDistance(approx, inf), ULPDistance(approx, sup));

		if (dist <= 5)
			continue;

		if (x > worst)
		{
			worst = x;
			std::cout << std::format("{}\n", worst);
		}
	}
	*/

	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ -0.2, -1e-300, false };

	MaxULPRounded(ReferenceLambertWM1, MakeSerial<MuirpairWm1>, [&]() { return dist(gen); });
}
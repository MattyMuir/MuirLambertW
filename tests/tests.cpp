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

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ EM_UP, -0.25 };

	double worst = -1;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = ReferenceLambertW0(x);

		double approx = MakeSerial<MuirpairW0>(x);

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
}
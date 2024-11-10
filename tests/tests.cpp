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
#include "ReferenceW.h"

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

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ EM_UP, -0.2 };

	double largestFailing = -INFINITY;
	for (;;)
	{
		double x = dist(gen);

		double approx = MakeSerial<MuirWm1>(x);
		auto exact = ReferenceWm1(x);
		if (ULPDistance(approx, exact) >= 5 && x > largestFailing)
		{
			std::cout << std::format("{}\n", x);
			largestFailing = x;
			dist = std::uniform_real_distribution<double>{ largestFailing, -0.2 };
		}
	}
}
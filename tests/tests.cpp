#include <iostream>
#include <vector>
#include <random>
#include <format>
#include <chrono>

#include <MuirLambertW.h>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <boost/math/special_functions/lambert_w.hpp>

#include "constants.h"

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

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };

#if 1
	ReciprocalDistributionEx<double> dist{ EM_UP, 0, false };
	MaxULPRounded(ReferenceW0, [](double x) { return MuirW0(x); }, [&]() { return dist(gen); }, 0);
#else
	std::uniform_real_distribution<double> dist{ EM_UP, -0.2 };

	for (;;)
	{
		double x = dist(gen);

		auto approx = MuirW0(x);
		auto exact = ReferenceW0(x);

		if (ULPDistance(approx, exact) > 4)
		{
			dist = std::uniform_real_distribution<double>{ x, -0.2 };
			std::cout << std::format("{}\n", x);
		}
	}
#endif
}
#include <iostream>
#include <vector>
#include <random>
#include <format>
#include <chrono>

#include <MuirLambertW.h>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <boost/math/special_functions/lambert_w.hpp>

#include "../bench/others/BarryLambertW.h"
#include "../bench/others/VebericLambertW.h"

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

double ExpMapW0(double x)
{
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	if (x < 700)
		return EM_UP / (1 + exp(x));

	return EM_UP / (1 + exp(x - 62)) / 8.4383566687414544891e+26;
}

float ExpMapWm1(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP / (1 + exp(x));
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<double> dist{ EM_UP, INFINITY, false };

	//MaxULPRounded(ReferenceW0, [](double x) { return MuirW0v2(x); }, []() { return dist(gen); }, 0);
	freopen("err.csv", "w", stdout);
	ULPHistogram(ReferenceW0, [](double x) { return MuirW0v2(x); }, -0.365, -0.05, 0.005, IdentityMap, 100'000);
}
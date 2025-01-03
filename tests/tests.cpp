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
#include "../bench/others/FukushimaLambertW.h"

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

float ExpMapW0(float x)
{
	return EM_UPf + exp(x);
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

template <typename Ty>
using MapTy = Ty(*)(Ty);

template <typename Ty>
std::pair<std::vector<Ty>, std::vector<UIntType<Ty>>> ULPHistogramVals(auto referenceFunc, auto approxFunc, Ty min, Ty max, Ty step, MapTy<Ty> map = IdentityMap, size_t iter = 10'000)
{
	static std::mt19937_64 gen{ std::random_device{}() };

	std::vector<Ty> xs;
	std::vector<UIntType<Ty>> errs;
	for (Ty low = min; low < max; low += step)
	{
		Ty high = low + step;
		std::uniform_real_distribution<Ty> dist{ low, high };

		UIntType<Ty> err = MaxULPRounded(referenceFunc, approxFunc, [&]() { return map(dist(gen)); }, iter);
		xs.push_back(std::midpoint(low, high));
		errs.push_back(err);
	}

	return { xs, errs };
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx<float> dist{ -0.051, 0.051, false };

	MaxULPRounded(ReferenceW0f, [](float x) { return MuirW0v2(x); }, [&]() { return dist(gen); }, 0);
}
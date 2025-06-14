#include <iostream>
#include <format>
#include <fstream>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <MuirLambertW.h>

#include "../bench/others/BarryLambertW.h"
#include "../bench/others/VebericLambertW.h"
#include "../bench/others/FukushimaLambertW.h"
#include "../bench/others/MuirFukushima.h"
#include "../bench/others/PsemLambertW.h"
#include "../bench/others/FukushimaMinimax.h"
#include <boost/math/special_functions/lambert_w.hpp>

#define TIMER_NPRINT
#include "../bench/Timer.h"

static constexpr float EM_UPf = -0.36787942f;
static constexpr double EM_UP = -0.3678794411714423;

Intervalf ReferenceW0f(float x) { thread_local ReferenceWf evaluator; return evaluator.W0(x); }
Intervalf ReferenceWm1f(float x) { thread_local ReferenceWf evaluator; return evaluator.Wm1(x); }
Interval ReferenceW0(double x) { thread_local ReferenceW evaluator; return evaluator.W0(x); }
Interval ReferenceWm1(double x) { thread_local ReferenceW evaluator; return evaluator.Wm1(x); }

float ExpMapW0(float x)
{
	return EM_UPf + exp(x);
}

float ExpMapWm1(float x)
{
	if (x < 88.7f)
		return EM_UPf / (1 + expf(x));

	return EM_UPf / (1 + exp(x - 14)) * 8.315287e-07;
}

double ExpMapW0(double x)
{
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	if (x < 700)
		return EM_UP / (1 + exp(x));

	return EM_UP / (1 + exp(x - 62)) * 1.185064864233981e-27;
}

using FloatTy = float;
std::vector<UIntType<FloatTy>> GetMaximumError(auto referenceFunc, auto approxFunc, const std::vector<FloatTy>& mins, const std::vector<FloatTy>& maxs)
{
	static constexpr size_t Iter = 10'000;
	static std::mt19937_64 gen{ std::random_device{}() };

	std::vector<UIntType<FloatTy>> errs;
	for (size_t binIdx = 0; binIdx < mins.size(); binIdx++)
	{
		UIntType<FloatTy> maxError = 0;
		std::uniform_real_distribution<FloatTy> dist{ mins[binIdx], maxs[binIdx] };
		for (size_t i = 0; i < Iter; i++)
		{
			FloatTy x = ExpMapWm1(dist(gen));

			auto exact = referenceFunc(x);
			FloatTy approx = approxFunc(x);
			UIntType<FloatTy> err = ULPDistance(approx, exact);
			maxError = std::max(maxError, err);
		}

		errs.push_back(maxError);
	}

	return errs;
}

static inline double AddEm(double x)
{
	static constexpr double emHigh = 0.36787944117144232160;
	static constexpr double emLow = -1.2428753672788363168e-17;

	return (x + emHigh) + emLow;
}

static inline float AddEmf(float x)
{
	return (x + 0.36787945f) - 9.149756e-09f;
}

double Approx(double x)
{
	static constexpr double P[] = {
		-120289.93485476710657323,3.327077642475286357704e6,-1.9206164555096615068559e6,518928.9216152833057284,832465.5446416577812005,-308960.153561595343874,15325.495522621048313411,41424.65862966648660759,4678.698788664447341361,140.69089035455689735618,1.0000080595936066116505
	};
	static constexpr double Q[] = {
		0, 97199.08977186529421392,332447.50807614386759464,-91215.91176509049906289,90190.20689480083183679,47433.780149220537699849,4832.1871093625870458303,141.6974875853895953604,1.
	};

	double t = sqrt(AddEm(x));

	double numer = P[9];
	for (size_t i = 0; i < 9; i++)
		numer = numer * t + P[8 - i];

	double denom = Q[8];
	for (size_t i = 0; i < 8; i++)
		denom = denom * t + Q[7 - i];

	return x / (numer / denom + 0.375 + t * 1.5);
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<double> dist{ EM_UP, INFINITY, false};

	ErrorSearcher searcher{ ReferenceW0, Overload<double, MuirW0> };
	searcher.MaxError([]() { return dist(gen); });
}
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
		0.367879440562318100277058,0.01629579856329787014053,-1.95895520261828626809831,-2.41320385931197971594225,-0.86154030099045003352524,0.008792042127552119907556,0.0459800172381655582073877,0.00552259547416009664176818,0.00014590301513459501193605,2.536810387426921991311715e-8
	};
	static constexpr double Q[] = {
		1.,3.149216216684359316991575,3.546977408786420079754295,1.775881710680574587670422,0.4031482003338336551181654,0.03796487558820609631875625,0.001139529814620952218770483,2.909604174184339899208132e-6,-1.268405193713460995655858e-8
	};

	double t = sqrt(AddEm(x));

	double numer = P[9];
	for (size_t i = 0; i < 9; i++)
		numer = numer * t + P[8 - i];

	double denom = Q[8];
	for (size_t i = 0; i < 8; i++)
		denom = denom * t + Q[7 - i];

	return x / (t * 2 + numer / denom);
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<double> dist{ EM_UP, 0, false};

	ErrorSearcher searcher{ ReferenceWm1, MuirWm1v2 };
	searcher.MaxError([]() { return dist(gen); });
}
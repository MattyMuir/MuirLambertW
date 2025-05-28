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

Intervalf ReferenceW0f(float x) { static ReferenceWf evaluator; return evaluator.W0(x); }
Intervalf ReferenceWm1f(float x) { static ReferenceWf evaluator; return evaluator.Wm1(x); }
Interval ReferenceW0(double x) { static ReferenceW evaluator; return evaluator.W0(x); }
Interval ReferenceWm1(double x) { static ReferenceW evaluator; return evaluator.Wm1(x); }

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

float AddEmf(float x)
{
	return (x + 0.36787945f) - 9.149756e-09f;
}

float Approx(float x)
{
	static constexpr double e2 = 5.43656365691809;
	double p = sqrt(e2 * x + 2.0);

	static constexpr double P[] = {
		-1.0000000289164983219173,0.7484377282132577749822601,0.3718673404409796721564797,-0.3760756806217422409838732,0.06188486908948103343472952
	};

	static constexpr double Q[] = {
		1,-1.748435627144009813001159,1.043209266905849599870069,-0.2369809021813308821994934,0.0145973490111531148723476
	};

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * p + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * p + Q[3 - i];

	return numer / denom;
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<float> dist{ EM_UPf, INFINITY, false };

	MaxULPRounded(ReferenceW0f, MuirW0v2, []() { return dist(gen); }, 0);
}
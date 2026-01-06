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

float Approx(float x)
{
	static constexpr float P[] = {
		2103.538895966433342674,-18492.4787689402541616248,4428.6774419532113647695,-171.385458185158872416365,1.00158289237695144114984
	};
	static constexpr float Q[] = {
		-9738.4370098199988067853,3751.23599762249746210982,-165.078029162292346655659,1.
	};

	float t = log(-x);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * t + Q[2 - i];

	return numer / denom;
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<float> dist{ EM_UPf, 0.0f, false};

	ErrorSearcher searcher{ ReferenceWm1f, MakeSerial<float, MuirWm1> };
	searcher.MaxError([]() { return dist(gen); });
}
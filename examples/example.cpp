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
		2.527351077301677171985851,1.376271973918832998774399,0.061118686508232127572727,-0.141295654553723696451935,-0.033728536316733169831178,-0.0025420185447289137885008,-0.00006664589607946860036686,-4.95157762697489553558e-7,-2.5731103350149386285e-10,2.587497211686779031e-13
	};
	static constexpr double Q[] = {
		1.,1.104130249984394722360317,0.5100993866359004287957678,0.1150243579777564561094546,0.01216796645784279085341347,0.0005632651755312496245963974,0.00001005684306251249387739756,5.02779367612422010297913e-8
	};

	double t = sqrt(-1 - log(-x));

	double numer = P[9];
	for (size_t i = 0; i < 9; i++)
		numer = numer * t + P[8 - i];

	double denom = Q[7];
	for (size_t i = 0; i < 7; i++)
		denom = denom * t + Q[6 - i];

	return (numer / denom - 3.52734375) - t * t;
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<double> dist{ -0.0000453999297625, 0, false};

	ErrorSearcher searcher{ ReferenceWm1, Approx };
	searcher.MaxError([]() { return dist(gen); });
}
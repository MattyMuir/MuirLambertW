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
	double t = x / (1 + x);

	static constexpr double P[] = {
		3.821005607016648530645815,-0.835427227384345479175442,-0.1676110305829701516749124,0.0000697878890147964021931,0.0002706420731114464410675,4.073286490834762983453e-6,1.2696378656088000085754e-8,2.594658242216529757192e-12,-8.34652818821295159020687e-16
	};

	static constexpr double Q[] = {
		1.,0.01279854337665376068872497,-0.0964334497380958423952366,-0.008782032545023475891919091,-0.0002154194241684539978467438,-1.62242239699663351508148e-6,-2.99518766936108213590929e-9,8.346528188212951590206866e-16
	};

	double numer = P[6];
	for (size_t i = 0; i < 6; i++)
		numer = numer * t + P[5 - i];

	double denom = Q[6];
	for (size_t i = 0; i < 6; i++)
		denom = denom * t + Q[5 - i];

	return t * (1 + numer / denom);
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<double> dist{ EM_UP, INFINITY, false};

	ErrorSearcher searcher{ ReferenceW0, MuirW0v2 };
	searcher.MaxError([]() { return dist(gen); });
}
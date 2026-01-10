#include <iostream>
#include <format>
#include <fstream>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <MuirLambertW.h>

#include "../bench/others/BarryLambertW.h"
#include "../bench/others/VebericLambertWOld.h"
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
		-34120.78392669408956,5199.530423442335303,63.39631437058543741,-2.5877024603665257607,0.0015828923769514411498
	};
	static constexpr float Q[] = {
		-9738.437009819998806785,3751.2359976224974621098,-165.07802916229234665566,1.
	};

	float t = logf(-x);

	float numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	float denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * t + Q[2 - i];

	return numer / denom + t - 3.7197265625f;
}

template <typename Ty>
using MapTy = Ty(*)(Ty);

template <typename Ty>
std::vector<UIntType<Ty>> ULPHistogramVals(auto referenceFunc, auto approxFunc, Ty min, Ty max, Ty step, MapTy<Ty> map = IdentityMap, size_t iter = 10'000)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::vector<UIntType<Ty>> errs;

	for (Ty low = min; low < max; low += step)
	{
		Ty high = low + step;
		std::uniform_real_distribution<Ty> dist{ low, high };

		UIntType<Ty> err = MaxULPRounded(referenceFunc, approxFunc, [&]() { return map(dist(gen)); }, iter);
		errs.push_back(err);
	}

	return errs;
}

int main()
{
	std::vector<uint64_t> proposed = ULPHistogramVals(ReferenceW0, Overload<double, MuirW0>, -40.0, 0.0, 1.0, ExpMapW0);

	std::vector<std::vector<uint64_t>> others;
	others.push_back(ULPHistogramVals(ReferenceW0, BarryLambertW0,                           -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, veberic_old::LambertW<0>,                 -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, utl::LambertW<0>,                         -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, boost::math::lambert_w0<double>,          -40.0, 0.0, 1.0, ExpMapW0));
	//others.push_back(ULPHistogramVals(ReferenceW0, PsemLambertW0,                            -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, Fukushima::LambertW0,                     -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, Overload<double, MuirFukushimaW0>,        -40.0, 0.0, 1.0, ExpMapW0));
	others.push_back(ULPHistogramVals(ReferenceW0, Overload<double, FukushimaMinimaxW0>,     -40.0, 0.0, 1.0, ExpMapW0));

	std::vector<uint64_t> mins = others[0];
	size_t numRows = others[0].size();
	size_t numCols = others.size();
	for (size_t r = 0; r < numRows; r++)
		for (size_t c = 1; c < numCols; c++)
			mins[r] = std::min(mins[r], others[c][r]);

	for (size_t r = 0; r < numRows; r++)
		std::cout << std::format("{:.2f} {} {}\n", -40.0 + 1.0 * r, mins[r], proposed[r]);
}
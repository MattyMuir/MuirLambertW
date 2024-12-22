#include <iostream>
#include <vector>
#include <format>
#include <random>
#include <fstream>

#include <immintrin.h>

#define TIMER_NPRINT
#include "Timer.h"

// Implementations
#include <MuirLambertW.h>
#include "others/FukushimaLambertW.h"
#include "others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"
#include "others/VebericLambertW.h"

using Function1D = double(*)(double);
using SimdFunction1D = __m256d(*)(__m256d);

void ApplyFunction(Function1D func, std::vector<double>& dst, const std::vector<double>& src)
{
	for (size_t i = 0; i < src.size(); i++)
		dst[i] = func(src[i]);
}

void ApplyFunction(SimdFunction1D func, std::vector<double>& dst, const std::vector<double>& src)
{
	for (size_t i = 0; i < src.size(); i += 4)
	{
		__m256d vals = _mm256_loadu_pd(src.data() + i);
		__m256d res = func(vals);
		_mm256_storeu_pd(dst.data() + i, res);
	}
}

double TimeFunction(const auto& function, const std::vector<double>& src)
{
	std::vector<double> dst(src.size());

	Timer t;
	ApplyFunction(function, dst, src);
	t.Stop();

	return t.GetSeconds();
}

std::vector<double> CreateArray(size_t size, double min, double max)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ min, max };

	std::vector<double> ret;
	ret.reserve(size);

	for (size_t i = 0; i < size; i++)
		ret.push_back(dist(gen));

	return ret;
}

double ExpMapW0(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	static constexpr double EM_UP = -0.3678794411714423;
	return EM_UP / (1 + exp(x));
}

int main()
{
	// === Parameters ===
	static constexpr size_t ArrSize = 1'200;
	static constexpr size_t Repeats = 420;
	double binMin = -3;
	double binMax = 20;
	double binWidth = 0.2;
	// ==================

	std::ofstream file{ "arraybench.csv" };

	file << "Min,Max,Barry,Veberic,Fukushima,Boost,Muir,MuirSerial\n";
	for (double min = binMin; min < binMax; min += binWidth)
	{
		double max = min + binWidth;
		std::vector<double> src = CreateArray(ArrSize, ExpMapW0(min), ExpMapW0(max));

		double barryTime = 0, vebericTime = 0, fukushimaTime = 0, boostTime = 0, muirTime = 0, muirSerialTime = 0;
		for (size_t repeat = 0; repeat < Repeats; repeat++)
		{
			//barryTime += TimeFunction(BarryLambertW0, src);
			//vebericTime += TimeFunction(utl::LambertW<0>, src);
			fukushimaTime += TimeFunction(Fukushima::LambertW0, src);
			//boostTime += TimeFunction(boost::math::lambert_w0<double>, src);
			//muirTime += TimeFunction([](__m256d x) { return MuirW0(x); }, src);
			//muirSerialTime += TimeFunction([](double x) { return MuirW0(x); }, src);
		}

		barryTime /= Repeats;
		vebericTime /= Repeats;
		fukushimaTime /= Repeats;
		boostTime /= Repeats;
		muirTime /= Repeats;
		muirSerialTime /= Repeats;

		file << std::format("{:.2f},{:.2f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f}\n", min, max, barryTime, vebericTime, fukushimaTime, boostTime, muirTime, muirSerialTime);
		std::cout << min << " - " << max << '\n';
	}
}
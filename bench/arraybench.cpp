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
#include <boost/math/special_functions/lambert_w.hpp>
#include "others/FukushimaLambertW.h"
#include "others/BarryLambertW.h"
#include "others/VebericLambertW.h"
#include "others/VebericLambertWOld.h"
#include "others/MuirFukushima.h"
#include "others/PsemLambertW.h"
#include "others/FukushimaMinimax.h"

#define FORMAT_CSV 0

using Function1D = double(*)(double);
using SimdFunction1D = __m256d(*)(__m256d);

static constexpr double EM_UP = -0.3678794411714423;

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
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	if (x < 700)
		return EM_UP / (1 + exp(x));

	return EM_UP / (1 + exp(x - 62)) * 1.185064864233981e-27;
}

int main()
{
	// === Parameters ===
	static constexpr size_t ArrSize = 256;
	static constexpr size_t Repeats = 10'000;
	double binMin = -30;
	double binMax = 30;
	size_t binNum = 500;
	size_t benchNum = 10;
	double binWidth = (binMax - binMin) / binNum;
	bool UseThroughput = false;
	// ==================

#if FORMAT_CSV
	std::ofstream file{ "arraybench.csv" };
	file << "min,max,barry,veberic,vebericold,fukushima,boost,muir,muirserial,muirfukushima,psem,fukushimaminimax\n";
#else
	std::ofstream file{ "arraybench.dat" };
	file << "min max barry veberic vebericold fukushima boost muir muirserial muirfukushima psem fukushimaminimax\n";
#endif

	std::vector<std::vector<double>> timings(binNum, std::vector<double>(benchNum));

	for (size_t repeat = 0; repeat < Repeats; repeat++)
	{
		for (size_t binIdx = 0; binIdx < binNum; binIdx++)
		{
			// Get reference to timing bin
			std::vector<double>& binTimings = timings[binIdx];

			// Create array
			double min = binMin + binIdx * binWidth;
			double max = binMin + (binIdx + 1) * binWidth;
			std::vector<double> src = CreateArray(ArrSize, ExpMapWm1(min), ExpMapWm1(max));

			// Time functions
			binTimings[0] += TimeFunction(BarryLambertWm1, src);
			binTimings[1] += TimeFunction(utl::LambertW<-1>, src);
			binTimings[2] += TimeFunction(veberic_old::LambertW<-1>, src);
			binTimings[3] += TimeFunction(Fukushima::LambertWm1, src);
			binTimings[4] += TimeFunction(boost::math::lambert_wm1<double>, src);
			binTimings[5] += TimeFunction([](__m256d x) { return MuirWm1(x); }, src);
			binTimings[6] += TimeFunction([](double x) { return MuirWm1(x); }, src);
			binTimings[7] += TimeFunction([](double x) { return MuirFukushimaWm1(x); }, src);
			binTimings[8] += TimeFunction(PsemLambertWm1, src);
			binTimings[9] += TimeFunction([](double x) { return FukushimaMinimaxWm1(x); }, src);
		}

		std::cout << repeat << '\n';
	}

	char sep = FORMAT_CSV ? ',' : ' ';
	for (size_t binIdx = 0; binIdx < binNum; binIdx++)
	{
		// Get reference to timing bin
		std::vector<double>& binTimings = timings[binIdx];

		// Print bin range
		double min = binMin + binIdx * binWidth;
		double max = binMin + (binIdx + 1) * binWidth;
		file << std::format("{:.4}{}{:.4}", min, sep, max);

		for (double time : binTimings)
		{
			if (UseThroughput)
				file << std::format("{}{:.4}", sep, Repeats / time * ArrSize * 1e-6);
			else
				file << std::format("{}{:.10}", sep, time / Repeats);
		}
			
		file << '\n';
	}
}
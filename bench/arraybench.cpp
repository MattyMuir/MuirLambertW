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

// === Parameters ===
static constexpr size_t ArrSize = 256;
static constexpr size_t Repeats = 100;
double binMin = -30;
double binMax = 30;
size_t binNum = 500;
size_t benchNum = 10;
double binWidth = (binMax - binMin) / binNum;
bool UseThroughput = false;
// ==================

#define RESET TimeFunction([](double x) { return sqrt(x); }, src)

#define MAKE_TEMPLATED_WRAPPER(name, func, targs)\
template <int64_t Branch>\
double name##_wrapper(double x)\
{\
	if constexpr (Branch == 0)\
		return func##0 targs(x);\
	if constexpr (Branch == -1)\
		return func##m1 targs(x);\
}

#define MAKE_TEMPLATED_WRAPPER_SIMD(name, func, targs)\
template <int64_t Branch>\
__m256d name##_wrapper(__m256d x)\
{\
	if constexpr (Branch == 0)\
		return func##0 targs(x);\
	if constexpr (Branch == -1)\
		return func##m1 targs(x);\
}

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

void WriteToFile(const std::string& filepath, const std::vector<std::vector<double>>& timings)
{
	std::cout << std::format("Saving results to {}\n", filepath);
	std::ofstream file{ filepath };
	file << "min max barry veberic vebericold fukushima boost muir muirserial muirfukushima psem fukushimaminimax\n";

	for (size_t binIdx = 0; binIdx < binNum; binIdx++)
	{
		// Print bin range
		const std::vector<double>& binTimings = timings[binIdx];
		float min = binMin + binIdx * binWidth;
		file << std::format("{:.4}{}{:.4}", min, ' ', min + binWidth);

		for (double time : binTimings)
		{
			double result = UseThroughput ? (Repeats / time * ArrSize * 1e-6) : (time / Repeats);
			file << std::format("{}{:.10}", ' ', result);
		}
		file << '\n';
	}
}

MAKE_TEMPLATED_WRAPPER(Barry, BarryLambertW, )
MAKE_TEMPLATED_WRAPPER(Fukushima, Fukushima::LambertW, )
MAKE_TEMPLATED_WRAPPER(Boost, boost::math::lambert_w, <double>)
MAKE_TEMPLATED_WRAPPER(Psem, PsemLambertW, )
MAKE_TEMPLATED_WRAPPER(MuirFukushima, MuirFukushimaW, )
MAKE_TEMPLATED_WRAPPER(FukushimaMinimax, FukushimaMinimaxW, )
MAKE_TEMPLATED_WRAPPER(ExpMap, ExpMapW, )
MAKE_TEMPLATED_WRAPPER(Muir, MuirW, )
MAKE_TEMPLATED_WRAPPER_SIMD(MuirSimd, MuirW, )

template <int64_t Branch>
std::vector<std::vector<double>> RunBenchmark()
{
	std::vector<std::vector<double>> timings(binNum, std::vector<double>(benchNum));

	for (size_t repeat = 0; repeat < Repeats; repeat++)
	{
		if (repeat % 10 == 0)
			std::cout << std::format("Progress: {:.3f}%\r", (double)repeat / Repeats * 100.0);

		for (size_t binIdx = 0; binIdx < binNum; binIdx++)
		{
			// Get reference to timing bin
			std::vector<double>& binTimings = timings[binIdx];

			// Create array
			double min = binMin + binIdx * binWidth;
			double max = binMin + (binIdx + 1) * binWidth;
			std::vector<double> src = CreateArray(ArrSize, ExpMap_wrapper<Branch>(min), ExpMap_wrapper<Branch>(max));

			// Time functions
			RESET;
			binTimings[0] += TimeFunction(Barry_wrapper<Branch>, src);+
			RESET;
			binTimings[1] += TimeFunction(utl::LambertW<Branch>, src);
			RESET;
			binTimings[2] += TimeFunction(veberic_old::LambertW<Branch>, src);
			RESET;
			binTimings[3] += TimeFunction(Fukushima_wrapper<Branch>, src);
			RESET;
			binTimings[4] += TimeFunction(Boost_wrapper<Branch>, src);
			RESET;
			binTimings[5] += TimeFunction(MuirSimd_wrapper<Branch>, src);
			RESET;
			binTimings[6] += TimeFunction(Muir_wrapper<Branch>, src);
			RESET;
			binTimings[7] += TimeFunction(MuirFukushima_wrapper<Branch>, src);
			RESET;
			binTimings[8] += TimeFunction(Psem_wrapper<Branch>, src);
			RESET;
			binTimings[9] += TimeFunction(FukushimaMinimax_wrapper<Branch>, src);
			RESET;
		}
	}

	std::cout << "Finished!         \n";
	return timings;
}

int main(int argc, char** argv)
{
	if (argc != 2)
	{
		std::cout << "Provide an argument specifying which branch to benchmark (0/-1)\n";
		return -1;
	}

	int64_t branch = std::stoll(argv[1]);

	std::vector<std::vector<double>> timings;
	if (branch == 0) timings = RunBenchmark<0>();
	else if (branch == -1) timings = RunBenchmark<-1>();
	else { std::cout << "Invalid branch!\n"; return -1; }

	WriteToFile(std::format("arraybench{}.dat", branch ? "m1" : "0"), timings);
}
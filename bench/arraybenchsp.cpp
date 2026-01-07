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
static constexpr size_t Repeats = 2'000;
float binMin = -15.0f;
float binMax = 15.0f;
size_t binNum = 500;
size_t benchNum = 10;
float binWidth = (binMax - binMin) / binNum;
bool UseThroughput = false;
// ==================

#define RESET TimeFunction(sqrtf, src)

#define MAKE_TEMPLATED_WRAPPER(name, func, targs)\
template <int64_t Branch>\
float name##_wrapper(float x)\
{\
	if constexpr (Branch == 0)\
		return func##0 targs(x);\
	if constexpr (Branch == -1)\
		return func##m1 targs(x);\
}

#define MAKE_TEMPLATED_WRAPPER_SIMD(name, func, targs)\
template <int64_t Branch>\
__m256 name##_wrapper(__m256 x)\
{\
	if constexpr (Branch == 0)\
		return func##0 targs(x);\
	if constexpr (Branch == -1)\
		return func##m1 targs(x);\
}

using Function1Df = float(*)(float);
using SimdFunction1Df = __m256(*)(__m256);

static constexpr float EM_UPf = -0.36787942f;

void ApplyFunction(Function1Df func, std::vector<float>& dst, const std::vector<float>& src)
{
	for (size_t i = 0; i < src.size(); i++)
		dst[i] = func(src[i]);
}

void ApplyFunction(SimdFunction1Df func, std::vector<float>& dst, const std::vector<float>& src)
{
	for (size_t i = 0; i < src.size(); i += 8)
	{
		__m256 vals = _mm256_loadu_ps(src.data() + i);
		__m256 res = func(vals);
		_mm256_storeu_ps(dst.data() + i, res);
	}
}

double TimeFunction(const auto& function, const std::vector<float>& src)
{
	std::vector<float> dst(src.size());

	Timer t;
	ApplyFunction(function, dst, src);
	t.Stop();

	return t.GetSeconds();
}

std::vector<float> CreateArray(size_t size, float min, float max)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<float> dist{ min, max };

	std::vector<float> ret;
	ret.reserve(size);

	for (size_t i = 0; i < size; i++)
		ret.push_back(dist(gen));

	return ret;
}

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

float VebericFloatW0(float x) { return utl::LambertW<0>(x); }
float VebericFloatWm1(float x) { return utl::LambertW<-1>(x); }
float VebericOldFloatW0(float x) { return veberic_old::LambertW<0>(x); }
float VebericOldFloatWm1(float x) { return veberic_old::LambertW<-1>(x); }

MAKE_TEMPLATED_WRAPPER(Barry, BarryLambertW,)
MAKE_TEMPLATED_WRAPPER(Veberic, VebericFloatW,)
MAKE_TEMPLATED_WRAPPER(VebericOld, VebericOldFloatW,)
MAKE_TEMPLATED_WRAPPER(Fukushima, Fukushima::LambertW,)
MAKE_TEMPLATED_WRAPPER(Boost, boost::math::lambert_w, <float>)
MAKE_TEMPLATED_WRAPPER(Psem, PsemLambertW,)
MAKE_TEMPLATED_WRAPPER(MuirFukushima, MuirFukushimaW,)
MAKE_TEMPLATED_WRAPPER(FukushimaMinimax, FukushimaMinimaxW,)
MAKE_TEMPLATED_WRAPPER(ExpMap, ExpMapW,)
MAKE_TEMPLATED_WRAPPER(Muir, MuirW,)
MAKE_TEMPLATED_WRAPPER_SIMD(MuirSimd, MuirW,)

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
			float min = binMin + binIdx * binWidth;
			float max = binMin + (binIdx + 1) * binWidth;
			std::vector<float> src = CreateArray(ArrSize, ExpMap_wrapper<Branch>(min), ExpMap_wrapper<Branch>(max));

			// Time functions
			RESET;
			binTimings[0] += TimeFunction(Barry_wrapper<Branch>, src);
			RESET;
			binTimings[1] += TimeFunction(Veberic_wrapper<Branch>, src);
			RESET;
			binTimings[2] += TimeFunction(VebericOld_wrapper<Branch>, src);
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

	WriteToFile(std::format("arraybenchf{}.dat", branch ? "m1" : "0"), timings);
}
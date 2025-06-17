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
#include "others/FukushimaMinimax.h"

#define FORMAT_CSV 1
#define RESET TimeFunction(sqrtf, src)

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

float BarryW0MadeFloat(float f)
{
	return BarryLambertW0(f);
}

float VebericW0MadeFloat(float f)
{
	return utl::LambertW<0>(f);
}

float VebericOldW0MadeFloat(float f)
{
	return veberic_old::LambertW<0>(f);
}

float BarryWm1MadeFloat(float f)
{
	return BarryLambertWm1(f);
}

float VebericWm1MadeFloat(float f)
{
	return utl::LambertW<-1>(f);
}

float VebericOldWm1MadeFloat(float f)
{
	return veberic_old::LambertW<-1>(f);
}

int main()
{
	// === Parameters ===
	static constexpr size_t ArrSize = 256;
	static constexpr size_t Repeats = 50'000;
	float binMin = -15.0f;
	float binMax = 15.0f;
	size_t binNum = 500;
	size_t benchNum = 10;
	float binWidth = (binMax - binMin) / binNum;
	bool UseThroughput = false;
	// ==================

#if FORMAT_CSV
	std::ofstream file{ "arraybench.csv" };
	file << "min,max,barry,veberic,vebericold,fukushima,boost,muir,muirserial,muirfukushima,fukushimaminimax,muirserialv2\n";
#else
	std::ofstream file{ "arraybench.dat" };
	file << "min max barry veberic vebericold fukushima boost muir muirserial muirfukushima fukushimaminimax muirserialv2\n";
#endif

	std::vector<std::vector<double>> timings(binNum, std::vector<double>(benchNum));

	for (size_t repeat = 0; repeat < Repeats; repeat++)
	{
		for (size_t binIdx = 0; binIdx < binNum; binIdx++)
		{
			// Get reference to timing bin
			std::vector<double>& binTimings = timings[binIdx];

			// Create array
			float min = binMin + binIdx * binWidth;
			float max = binMin + (binIdx + 1) * binWidth;
			std::vector<float> src = CreateArray(ArrSize, ExpMapW0(min), ExpMapW0(max));

			// Time functions
			//RESET;
			//binTimings[0] += TimeFunction(BarryW0MadeFloat, src);
			//RESET;
			//binTimings[1] += TimeFunction(VebericW0MadeFloat, src);
			//RESET;
			//binTimings[2] += TimeFunction(VebericOldW0MadeFloat, src);
			//RESET;
			//binTimings[3] += TimeFunction(Fukushima::LambertW0, src);
			RESET;
			binTimings[4] += TimeFunction(boost::math::lambert_w0<float>, src);
			RESET;
			binTimings[5] += TimeFunction([](__m256 x) { return MuirW0(x); }, src);
			RESET;
			binTimings[6] += TimeFunction([](float x) { return MuirW0(x); }, src);
			RESET;
			//binTimings[7] += TimeFunction([](float x) { return MuirFukushimaW0(x); }, src);
			//RESET;
			binTimings[8] += TimeFunction([](float x) { return FukushimaMinimaxW0(x); }, src);
			//RESET;
			//binTimings[9] += TimeFunction([](__m256 x) { return MuirW0(x); }, src);
			RESET;
		}

		std::cout << repeat << '\n';
	}

	char sep = FORMAT_CSV ? ',' : ' ';
	for (size_t binIdx = 0; binIdx < binNum; binIdx++)
	{
		// Get reference to timing bin
		std::vector<double>& binTimings = timings[binIdx];

		// Print bin range
		float min = binMin + binIdx * binWidth;
		float max = binMin + (binIdx + 1) * binWidth;
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
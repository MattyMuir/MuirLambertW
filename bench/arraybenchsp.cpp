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

using Function1D = float(*)(float);
using SimdFunction1D = __m256(*)(__m256);

void ApplyFunction(Function1D func, std::vector<float>& dst, const std::vector<float>& src)
{
	for (size_t i = 0; i < src.size(); i++)
		dst[i] = func(src[i]);
}

void ApplyFunction(SimdFunction1D func, std::vector<float>& dst, const std::vector<float>& src)
{
	for (size_t i = 0; i < src.size(); i += 8)
	{
		__m256 vals = _mm256_loadu_ps(src.data() + i);
		__m256 res = func(vals);
		_mm256_storeu_ps(dst.data() + i, res);
	}
}

float TimeFunction(const auto& function, const std::vector<float>& src)
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
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP + exp(x);
}

float ExpMapWm1(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP / (1 + exp(x));
}

int main()
{
	// === Parameters ===
	static constexpr size_t ArrSize = 1'000'000;
	static constexpr size_t Repeats = 30;
	float binMin = -10;
	float binMax = 10;
	float binWidth = 0.5;
	// ==================

	std::ofstream file{ "arraybench.csv" };

	file << "Min,Max,Fukushima,Barry,Boost,Muir\n";
	for (float min = binMin; min < binMax; min += binWidth)
	{
		float max = min + binWidth;
		std::vector<float> src = CreateArray(ArrSize, ExpMapWm1(min), ExpMapWm1(max));

		float fukushimaTime = 0, barryTime = 0, boostTime = 0, muirTime = 0;
		for (size_t repeat = 0; repeat < Repeats; repeat++)
		{
			//fukushimaTime += TimeFunction(Fukushima::LambertWm1, src);
			//barryTime += TimeFunction(BarryLambertWm1, src);
			boostTime += TimeFunction(boost::math::lambert_wm1<float>, src);
			muirTime += TimeFunction([](__m256 x) { return MuirWm1(x); }, src);
		}

		fukushimaTime /= Repeats;
		barryTime /= Repeats;
		boostTime /= Repeats;
		muirTime /= Repeats;

		file << std::format("{:.10f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f}\n", min, max, fukushimaTime, barryTime, boostTime, muirTime);
		std::cout << min << " - " << max << '\n';
	}
}
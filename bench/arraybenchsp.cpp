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

using Function1Df = float(*)(float);
using SimdFunction1Df = __m256(*)(__m256);

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
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP + exp(x);
}

float ExpMapWm1(float x)
{
	static constexpr float EM_UP = -0.36787942f;
	return EM_UP / (1 + exp(x));
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
	static constexpr size_t ArrSize = 1'000;
	static constexpr size_t Repeats = 100;
	float binMin = -10;
	float binMax = 20;
	float binWidth = 0.5;
	// ==================

	std::ofstream file{ "arraybench.csv" };

	file << "Min,Max,Barry,Veberic,VebericOld,Fukushima,Boost,Muir,MuirSerial\n";
	for (float min = binMin; min < binMax; min += binWidth)
	{
		float max = min + binWidth;
		std::vector<float> src = CreateArray(ArrSize, ExpMapW0(min), ExpMapW0(max));

		double barryTime = 0, vebericTime = 0, vebericOldTime = 0, fukushimaTime = 0, boostTime = 0, muirTime = 0, muirSerialTime = 0;
		for (size_t repeat = 0; repeat < Repeats; repeat++)
		{
			barryTime			+= TimeFunction(BarryW0MadeFloat, src);
			vebericTime			+= TimeFunction(VebericW0MadeFloat, src);
			vebericOldTime		+= TimeFunction(VebericOldW0MadeFloat, src);
			//fukushimaTime		+= TimeFunction(Fukushima::LambertW0, src);
			boostTime			+= TimeFunction(boost::math::lambert_w0<float>, src);
			muirTime			+= TimeFunction([](__m256 x) { return MuirW0(x); }, src);
			muirSerialTime		+= TimeFunction([](float x) { return MuirW0(x); }, src);
		}

		barryTime /= Repeats;
		vebericTime /= Repeats;
		vebericOldTime /= Repeats;
		fukushimaTime /= Repeats;
		boostTime /= Repeats;
		muirTime /= Repeats;
		muirSerialTime /= Repeats;

		file << std::format("{:.2f},{:.2f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f},{:.10f}\n", min, max, barryTime, vebericTime, vebericOldTime, fukushimaTime, boostTime, muirTime, muirSerialTime);
		std::cout << min << " - " << max << '\n';
	}
}
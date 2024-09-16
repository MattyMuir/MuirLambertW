#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"
#include "../tests/ReciprocalDistributionEx.h"

#include "others/FukushimaLambertW.h"
#include "others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"
#include "MuirLambertW.h"

#define BENCH_FUKUSHIMA 0
#define BENCH_BARRY 1
#define BENCH_BOOST 1
#define BENCH_MUIR_SIMD 1

#define BENCHMARK(func, name) _ += RunBenchmark(func, name)
#define SIMD_BENCHMARK(func, name) _2 = _mm256_add_pd(_2, RunBenchmark(func, name))

static constexpr double EM_UP = -0.3678794411714423;

// === Parameters ===
static constexpr size_t NumData = 1'000;
static constexpr size_t NumIter = 100'000;
// ==================

using BenchFunction = double(*)(double);
using SimdFunction = __m256d(*)(__m256d);

std::vector<double> data;
std::vector<__m256d> simdData;

void PrepareData()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx dist{ EM_UP, -std::numeric_limits<double>::min(), false};

	data.reserve(NumData);
	for (size_t i = 0; i < NumData; i++)
		data.push_back(dist(gen));

	simdData.reserve(NumData);
	for (size_t i = 0; i < NumData; i++)
		simdData.push_back(_mm256_set1_pd(dist(gen)));
}

double RunBenchmark(BenchFunction func, const char* name)
{
	double _ = 0.0;

	Timer t;
	for (size_t i = 0; i < NumIter; i++)
		for (size_t d = 0; d < NumData; d++)
			_ += func(data[d]);
	t.Stop();

	std::cout << std::format("{} took: {:.0f}ms\n", name, t.GetSeconds() * 1000);

	return _;
}

double MuirSimdMadeSerial(double x)
{
	return MuirWm1(_mm256_set1_pd(x))[0];
}

int main()
{
	PrepareData();

	double _ = 0.0;
	__m256d _2 = _mm256_setzero_pd();
#if BENCH_FUKUSHIMA
	BENCHMARK(Fukushima::LambertWm1, "Fukushima");
#endif
#if BENCH_BARRY
	BENCHMARK(BarryLambertWm1, "Barry");
#endif
#if BENCH_BOOST
	BENCHMARK(boost::math::lambert_wm1<double>, "Boost");
#endif
#if BENCH_MUIR_SIMD
	BENCHMARK(MuirSimdMadeSerial, "Muir SIMD");
#endif

	std::cout << _ << '\n';
	std::cout << _2[0];
}
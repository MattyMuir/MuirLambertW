#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"

#include "others/FukushimaLambertW.h"
#include "others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"
#include "MuirLambertW.h"

#define BENCH_FUKUSHIMA 0
#define BENCH_BARRY 0
#define BENCH_BOOST 1
#define BENCH_MUIR_SIMD 1

#define BENCHMARK(func, name) _ += RunBenchmark(func, name)

// === Parameters ===
static constexpr size_t NumData = 1'000;
static constexpr size_t NumIter = 100'000;
// ==================

using BenchFunction = double(*)(double);

std::vector<double> data;

void PrepareData()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static std::uniform_real_distribution<double> dist{ 1e44, 1e300 };

	data.reserve(NumData);
	for (size_t i = 0; i < NumData; i++)
		data.push_back(dist(gen));
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
	return MuirpairW0(_mm256_set1_pd(x))[0];
}

int main()
{
	PrepareData();

	double _ = 0.0;
#if BENCH_FUKUSHIMA
	BENCHMARK(Fukushima::LambertW0, "Fukushima");
#endif
#if BENCH_BARRY
	BENCHMARK(BarryLambertW0, "Barry");
#endif
#if BENCH_BOOST
	BENCHMARK(boost::math::lambert_w0<double>, "Boost");
#endif
#if BENCH_MUIR_SIMD
	BENCHMARK(MuirSimdMadeSerial, "Muir SIMD");
#endif

	std::cout << _;
}
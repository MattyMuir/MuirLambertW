#include <vector>
#include <random>

#define TIMER_NPRINT
#include "Timer.h"

#include "others/FukushimaLambertW.h"
#include "others/BarryLambertW.h"
#include "MuirLambertW.h"

#define BENCH_FUKUSHIMA 0
#define BENCH_BARRY 0
#define BENCH_MUIR 0
#define BENCH_MUIR_SIMD 1

void DoNotOptimize(auto... args)
{
	(std::cout << ... << args) << '\r';
	std::cout << "                                    \r";
}

int main()
{
	// === Parameters ===
	static constexpr size_t NumData = 1'000;
	static constexpr size_t NumIter = 100'000;
	// ==================

	// Prepare data
	static std::mt19937_64 gen{ std::random_device{}() };
	static std::uniform_real_distribution<double> dist{ 1, 1000 };

	std::vector<double> data;
	data.reserve(NumData);
	for (size_t i = 0; i < NumData; i++)
		data.push_back(dist(gen));

	double _0 = 0, _1 = 0, _2 = 0, _3 = 0;

#if BENCH_FUKUSHIMA
	TIMER(fukushima);
	for (size_t i = 0; i < NumIter; i++)
		for (size_t d = 0; d < NumData; d++)
			_0 += Fukushima::LambertW0(data[d]);
	STOP_LOG(fukushima);
#endif

#if BENCH_BARRY
	TIMER(barry);
	for (size_t i = 0; i < NumIter; i++)
		for (size_t d = 0; d < NumData; d++)
			_1 += BarryLambertW0(data[d]);
	STOP_LOG(barry);
#endif

#if BENCH_MUIR
	TIMER(muir);
	for (size_t i = 0; i < NumIter; i++)
		for (size_t d = 0; d < NumData; d++)
			_2 += MuirLambertW0(data[d]);
	STOP_LOG(muir);
#endif

#if BENCH_MUIR_SIMD
	TIMER(muir_simd);
	for (size_t i = 0; i < NumIter; i++)
		for (size_t d = 0; d < NumData; d++)
			_3 += MuirLambertW0Simd(_mm256_set1_pd(data[d]))[0];
	STOP_LOG(muir_simd);
#endif

	DoNotOptimize(_0, _1, _2, _3);
}
#include <iostream>
#include <format>

#include "MuirLambertW.h"

void PrintPack(__m256d x)
{
	double vals[4];
	_mm256_storeu_pd(vals, x);
	std::cout << std::format("[{}, {}, {}, {}]\n", vals[0], vals[1], vals[2], vals[3]);
}

int main()
{
	__m256d x = _mm256_setr_pd(3, 4, 5, 6);
	std::cout << "Input values: ";
	PrintPack(x);

	__m256d res = MuirW0(x);
	std::cout << "Results: ";
	PrintPack(res);
}
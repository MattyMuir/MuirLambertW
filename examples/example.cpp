#include <iostream>
#include <format>

#include "MuirLambertW.h"

int main()
{
	double res = MuirpairW0(_mm256_set1_pd(36))[0];
	std::cout << std::format("{}\n", res);
}
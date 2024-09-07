#include <iostream>

#include "MuirLambertW.h"

int main()
{
	std::cout << MuirLambertW0Simd(_mm256_set1_pd(3.0))[0];
}
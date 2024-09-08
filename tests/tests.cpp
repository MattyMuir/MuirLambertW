#include <random>

#include "ulp.h"
#include "wrappers.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"

int main(int argc, char** argv)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static std::uniform_real_distribution<double> dist{ 1, 1000 };

	MaxULPRounded(ReferenceLambertW0, MakeSerial<MuirLambertW0Simd>, [=]() { return dist(gen); });
}
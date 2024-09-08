#include <random>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"

int main(int argc, char** argv)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ EM_UP, -std::numeric_limits<double>::min(), false};
	
	MaxULPRounded(ReferenceLambertWM1, MakeSerial<MuirLambertWM1Simd>, [&]() { return dist(gen); });
}
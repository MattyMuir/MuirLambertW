#include <random>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"

double BoostLambertW0(double x)
{
	return boost::math::lambert_w0(x);
}

int main(int argc, char** argv)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ 0, Infinity, false};

	MaxULPRounded(ReferenceLambertW0, MakeSerial<MuirpairW0>, [&]() { return dist(gen); });
}
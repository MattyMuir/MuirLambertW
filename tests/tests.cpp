#include <random>
#include <format>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"

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
	// Last 5 (-0.33900 - -0.33800: 5) with fixed first coefficient
	// Last 5 (-0.34200 - -0.34100: 5)
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ EM_UP, Infinity, false };
	MaxULPRounded(ReferenceLambertW0, MakeSerial<MuirpairW0>, [&]() { return dist(gen); });

	//ULPHistogram(ReferenceLambertW0, MakeSerial<MuirpairW0>, -0.367, -0.330, 0.001);
}
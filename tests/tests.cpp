#include <vector>
#include <random>
#include <format>
#include <chrono>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ EM_UP, Infinity, false };

	MaxULPRounded(ReferenceLambertW0, MakeSerial<MuirW0>, [&]() { return dist(gen); });
}
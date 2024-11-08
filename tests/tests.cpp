#include <vector>
#include <random>
#include <format>
#include <chrono>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"
#include "constants.h"

#include "MuirLambertW.h"
#include "ReferenceW.h"

Interval ReferenceW0(double x)
{
	static ReferenceW evaluator;
	return evaluator.W0(x);
}

Interval ReferenceWm1(double x)
{
	static ReferenceW evaluator;
	return evaluator.Wm1(x);
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx dist{ EM_UP, -std::numeric_limits<double>::denorm_min(), false};

	// 14.776442463189579 has error of 5ULP
	MaxULPRounded(ReferenceWm1, MakeSerial<MuirWm1>, [&]() { return dist(gen); });
}
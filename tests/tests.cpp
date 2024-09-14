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

int main()
{
	ULPHistogram(ReferenceLambertW0, boost::math::lambert_w0, -20, 20, 0.5, [](double x) { return EM_UP + exp(x); });
}
#pragma once
#include "types.h"

double IdentityMap(double x);

using RandMap = double(*)(double);
void ULPHistogram(BoundedFunction1D boundedFunc, Function1D approxFunc, double min, double max, double step, RandMap map = IdentityMap);
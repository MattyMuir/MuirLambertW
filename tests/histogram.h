#pragma once
#include "types.h"

using RandMapf = float(*)(float);
float IdentityMap(float x);
void ULPHistogramSigned(BoundedFunction1Df boundedFunc, Function1Df approxFunc, float min, float max, float step, RandMapf map = IdentityMap);

using RandMap = double(*)(double);
double IdentityMap(double x);
void ULPHistogram(BoundedFunction1D boundedFunc, Function1D approxFunc, double min, double max, double step, RandMap map = IdentityMap);
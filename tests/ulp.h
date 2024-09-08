#pragma once
#include "types.h"

uint64_t ULPDistance(double a, double b);
uint64_t MaxULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, RandomFunction randFunc, uint64_t iter = 0);
double AvgULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, RandomFunction randFunc, uint64_t iter = 0);
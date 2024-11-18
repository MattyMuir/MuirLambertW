#pragma once
#include "types.h"

uint32_t ULPDistance(float a, float b);
int32_t ULPDistanceSigned(float a, float b);
uint32_t ULPDistance(float approx, Intervalf exact);
uint32_t MaxULPRounded(BoundedFunction1Df boundedFunc, Function1Df approxFunc, const RandomFunctionf& randFunc, uint64_t iter = 0);
std::pair<int32_t, int32_t> MaxULPSigned(BoundedFunction1Df boundedFunc, Function1Df approxFunc, const RandomFunctionf& randFunc, uint64_t iter = 0);

uint64_t ULPDistance(double a, double b);
uint64_t ULPDistance(double approx, Interval exact);
uint64_t MaxULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, const RandomFunction& randFunc, uint64_t iter = 0);
double AvgULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, const RandomFunction& randFunc, uint64_t iter = 0);
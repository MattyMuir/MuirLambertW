#include "ulp.h"

#include <iostream>
#include <format>
#include <bit>

#define LOGGING 1

uint64_t ULPDistance(double a, double b)
{
    if (a > b)
        std::swap(a, b);
    uint64_t aPunn = std::bit_cast<uint64_t>(a);
    uint64_t bPunn = std::bit_cast<uint64_t>(b);

    if (a >= 0 && b >= 0)
        return bPunn - aPunn;
    if (a < 0 && b < 0)
        return aPunn - bPunn;

    return (aPunn - 0x8000000000000000) + bPunn;
}

uint64_t MaxULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, RandomFunction randFunc, uint64_t iter)
{
    uint64_t maxError = 0;
	for (uint64_t i = 0; (i < iter) || !iter; i++)
	{
		double x = randFunc();

		double approx = approxFunc(x);
        if (!std::isfinite(approx))
            std::cout << std::format("NaN: {}\n", x);
		auto [inf, sup] = boundedFunc(x);

        uint64_t error = std::max(ULPDistance(approx, inf), ULPDistance(approx, sup));

        if (error > maxError)
        {
            maxError = error;
#if LOGGING
            std::cout << std::format("error: {} x: {}\n", maxError, x);
#endif
        }
	}

    return maxError;
}

double AvgULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, RandomFunction randFunc, uint64_t iter)
{
    double totalError = 0;
    for (uint64_t i = 0; (i < iter) || !iter; i++)
    {
        double x = randFunc();

        double approx = approxFunc(x);
        if (!std::isfinite(approx))
            continue;
        auto [inf, sup] = boundedFunc(x);

        uint64_t error = std::max(ULPDistance(approx, inf), ULPDistance(approx, sup));
        totalError += error;
    }

    return totalError / iter;
}
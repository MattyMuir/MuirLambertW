#include "ulp.h"

#include <iostream>
#include <format>
#include <bit>

#define LOGGING 0

uint32_t ULPDistance(float a, float b)
{
    if (a > b)
        std::swap(a, b);
    uint32_t aPunn = std::bit_cast<uint32_t>(a);
    uint32_t bPunn = std::bit_cast<uint32_t>(b);

    if (a >= 0 && b >= 0)
        return bPunn - aPunn;
    if (a < 0 && b < 0)
        return aPunn - bPunn;

    return (aPunn - 0x80000000) + bPunn;
}

int32_t ULPDistanceSigned(float a, float b)
{
    bool swapped = false;
    if (a > b)
    {
        std::swap(a, b);
        swapped = true;
    }
    uint32_t aPunn = std::bit_cast<uint32_t>(a);
    uint32_t bPunn = std::bit_cast<uint32_t>(b);

    uint32_t dist;
    if (a >= 0 && b >= 0)
        dist = bPunn - aPunn;
    else if (a < 0 && b < 0)
        dist = aPunn - bPunn;
    else
        dist = (aPunn - 0x80000000) + bPunn;

    return swapped ? -(int64_t)dist : (int64_t)dist;
}

uint32_t ULPDistance(float approx, Intervalf exact)
{
    return std::max(ULPDistance(approx, exact.inf), ULPDistance(approx, exact.sup));
}

uint32_t MaxULPRounded(BoundedFunction1Df boundedFunc, Function1Df approxFunc, const RandomFunctionf& randFunc, uint64_t iter)
{
    uint32_t maxError = 0;
    for (uint64_t i = 0; (i < iter) || !iter; i++)
    {
        // Sample random input
        float x = randFunc();

        // Evaluate approx and exact
        float approx = approxFunc(x);
        Intervalf exact = boundedFunc(x);

        if (!std::isfinite(approx) && std::isfinite(exact.inf) && std::isfinite(exact.sup))
            std::cout << std::format("Bad NaN: {}\n", x);

        uint32_t error = ULPDistance(approx, exact);

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

std::pair<int32_t, int32_t> MaxULPSigned(BoundedFunction1Df boundedFunc, Function1Df approxFunc, const RandomFunctionf& randFunc, uint64_t iter)
{
    int32_t maxErrorPos = INT32_MIN, minErrorNeg = INT32_MAX;
    for (uint64_t i = 0; (i < iter) || !iter; i++)
    {
        // Sample random input
        float x = randFunc();

        // Evaluate approx and exact
        float approx = approxFunc(x);
        Intervalf exact = boundedFunc(x);

        if (!std::isfinite(approx) && std::isfinite(exact.inf) && std::isfinite(exact.sup))
            std::cout << std::format("Bad NaN: {}\n", x);

        auto errInf = ULPDistanceSigned(approx, exact.inf);
        auto errSup = ULPDistanceSigned(approx, exact.sup);

        maxErrorPos = std::max(maxErrorPos, errInf);
        maxErrorPos = std::max(maxErrorPos, errSup);
        minErrorNeg = std::min(minErrorNeg, errInf);
        minErrorNeg = std::min(minErrorNeg, errSup);
    }

    return { minErrorNeg, maxErrorPos };
}

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

uint64_t ULPDistance(double approx, Interval exact)
{
    return std::max(ULPDistance(approx, exact.inf), ULPDistance(approx, exact.sup));
}

uint64_t MaxULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, const RandomFunction& randFunc, uint64_t iter)
{
    uint64_t maxError = 0;
	for (uint64_t i = 0; (i < iter) || !iter; i++)
	{
        // Sample random input
		double x = randFunc();

        // Evaluate approx and exact
		double approx = approxFunc(x);
        Interval exact = boundedFunc(x);

        if (!std::isfinite(approx) && std::isfinite(exact.inf) && std::isfinite(exact.sup))
            std::cout << std::format("Bad NaN: {}\n", x);

        uint64_t error = ULPDistance(approx, exact);

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

double AvgULPRounded(BoundedFunction1D boundedFunc, Function1D approxFunc, const RandomFunction& randFunc, uint64_t iter)
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
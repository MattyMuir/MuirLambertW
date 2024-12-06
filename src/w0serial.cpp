#include "w0serial.h"

#include <cstdint>
#include <cmath>
#include <bit>
#include <limits>

static inline double FirstApprox(double x)
{
    if (abs(x) < 1e-4)
        return x;

    static constexpr double P[] = {
        0,
        28.81244221654313,
        131.14615617307527,
        172.3834180355456,
        62.332230968272896,
        3.6088782677342017
    };

    static constexpr double Q[] = {
        28.813739177487342,
        159.9633664467659,
        289.069078387009,
        188.23555138648206,
        35.16238622673435,
        1
    };

    double numer = P[5];
    for (size_t i = 0; i < 5; i++)
        numer = numer * x + P[4 - i];

    double denom = Q[5];
    for (size_t i = 0; i < 5; i++)
        denom = denom * x + Q[4 - i];

    return numer / denom;
}

static inline double SecondApprox(double x)
{
    static constexpr double P[] = {
        64312.7454007891,
        43264.12227598657,
        20243.65384336377,
        453.17656235798086,
        1.0000432316050645
    };
    static constexpr double Q[] = {
        104342.57917932322,
        22499.368605590193,
        460.93750724715477,
        1
    };

    double logX = log(x);

    double numer = P[4];
    for (size_t i = 0; i < 4; i++)
        numer = numer * logX + P[3 - i];

    double denom = Q[3];
    for (size_t i = 0; i < 3; i++)
        denom = denom * logX + Q[2 - i];

    double approx = numer / denom;
    return approx;
}

static inline double GeneralW0(double x)
{
    double w = (x > 7.34) ? SecondApprox(x) : FirstApprox(x);

    // === Fritsch Iteration ===
    static constexpr double c23 = 2.0 / 3.0;
    double zn = log(x / w) - w;
    double temp = 1.0 + w;
    double temp2 = temp + c23 * zn;
    temp2 = 2.0 * temp * temp2;
    w = w * (1.0 + (zn / temp) * (temp2 - zn) / (temp2 - 2.0 * zn));

    return w;
}

static inline double AddEm(double x)
{
    static constexpr double emHigh = 0.36787944117144232160;
    static constexpr double emLow = -1.2428753672788363168e-17;

    return (x + emHigh) + emLow;
}

static inline double NearBranchSeries(double p)
{
    static constexpr double P[] = {
        -1.0,
        0.999999999999952,
        -0.3333333333272104,
        0.15277777746385043,
        -0.07962962100934082,
        0.04450217060974618,
        -0.02598313708983362,
        0.015623864373809608,
        -0.009555681655574338,
        0.0057905732001774755,
        -0.0032366081259528625,
        0.0014204580219827356,
        -0.00034378172909336254
    };

    // Evaluate polynomial using Horner's Method
    double value = P[12];
    for (size_t i = 0; i < 12; i++)
        value = value * p + P[11 - i];

    return value;
}

static inline double NearBranchW0(double x)
{
    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    return NearBranchSeries(p);
}

double MuirW0(double x)
{
    if (x == 0 || x == std::numeric_limits<double>::infinity())
        return x;

    return (x < -0.34) ? NearBranchW0(x) : GeneralW0(x);
}
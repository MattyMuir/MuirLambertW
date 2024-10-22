#include "w0serial.h"

#include <cstdint>
#include <cmath>
#include <bit>
#include <limits>

static inline double FirstApprox(double x)
{
    // === Constants ===
    static constexpr double e2 = 5.4365636569180904707;     // e * 2
    // =================

    if (abs(x) < 1e-4)
        return x;

    static constexpr double P[] = {
        -0.9998418255015216,
        0.9963519639934693,
        -0.31689383427598367,
        0.11927961836393368,
        -0.0389272342299362,
        0.009587149228046744,
        -0.0016619725895250465,
        0.000193387127759057,
        -1.426313511537818e-05,
        5.997840674383423e-07,
        -1.0923174871658654e-08
    };

    double reta = sqrt(x * e2 + 2);

    double approx = P[10];
    for (size_t i = 0; i < 10; i++)
        approx = approx * reta + P[9 - i];

    return approx;
}

static inline double SecondApprox(double x)
{
    static constexpr double P[] = {
        64393.137450661044568,
        43204.949550002405886,
        20295.724800471609342,
        453.37964270930132216,
        1.0
    };
    static constexpr double Q[] = {
        104344.40703457256313,
        22558.64516691800236,
        461.09954435682880103,
        0.9999372708768251572
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
    double w = (x > 20.0) ? SecondApprox(x) : FirstApprox(x);

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
        -1,
        0.9999999999999997,
        -0.3333333333330709,
        0.15277777774435128,
        -0.07962962801673715,
        0.044502274633980064,
        -0.025984120008670515,
        0.01562999136067536,
        -0.00958118607573627,
        0.00586072464198908,
        -0.003358979213132348,
        0.0015430644746481945,
        -0.0003975110776657685
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

    return (x < -0.3407) ? NearBranchW0(x) : GeneralW0(x);
}
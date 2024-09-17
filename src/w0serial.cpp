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
    // === Constants ===
    static constexpr double y = -5.70115661621093750e+00;
    static constexpr double offset = -0.132;
    // =================

    static constexpr double P[] = {
       6.42275660145116698e+00,
       1.33047964073367945e+00,
       6.72008923401652816e-02,
       1.16444069958125895e-03,
       7.06966760237470501e-06,
       5.48974896149039165e-09,
       -7.00379652018853621e-11,
       -1.89247635913659556e-13,
       -1.55898770790170598e-16,
    };
    static constexpr double Q[] = {
       1.00000000000000000e+00,
       3.34498588416632854e-01,
       2.51519862456384983e-02,
       6.81223810622416254e-04,
       7.94450897106903537e-06,
       4.30675039872881342e-08,
       1.10667669458467617e-10,
       1.31012240694192289e-13,
       6.53282047177727125e-17,
    };

    double logX = log(x + offset);

    double numer = P[8];
    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
    {
        numer = numer * logX + P[7 - i];
        denom = denom * logX + Q[7 - i];
    }

    double approx = numer / denom;
    approx = approx + y + logX;
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
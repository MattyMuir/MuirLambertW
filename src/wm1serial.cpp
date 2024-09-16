#include "wm1serial.h"

#include <cstdint>
#include <cmath>
#include <bit>
#include <limits>

#define GREATER 0x1E

double LogFast(double x)
{
    // === Constants ===
    static constexpr double ln2 = 0.69314718055994530942;
    static constexpr double dblMin = std::numeric_limits<double>::min();
    static constexpr double denormScale = 4503599627370496.0;               // 2^52
    static constexpr double denormOffset = 36.043653389117156090;           // ln(2^52)
    // =================

    // Fix for denormalized values
    bool isDenorm = x < dblMin;
    if (isDenorm)
        x *= denormScale;

    // Extract exponent and mantissa
    uint64_t punn = std::bit_cast<uint64_t>(x);
    int64_t exp = ((punn >> 52) - 1023);
    double mantissa = std::bit_cast<double>((punn & 0x3fffffffffffffff) | 0x3ff0000000000000);

    // Compute approximation
    static constexpr double P[] = {
        -1.7289291561920494e+00,
        2.78901155791566960e+00,
        -1.44093748876198707e+00,
        4.36015488686681152e-01,
        -5.50266844824230939e-02
    };

    double approx = P[4];
    for (size_t i = 0; i < 4; i++)
        approx = approx * mantissa + P[3 - i];

    // Fix for denormalized values
    if (isDenorm)
        approx -= denormOffset;

    return (double)exp * ln2 + approx;
}

double Approx(double x)
{
#if 0
    // === Constants ===
    static constexpr double s2 = 1.4142135623730950488;     // sqrt(2)
    static constexpr double c13 = 1.0 / 3.0;
    // =================

    double a =  (x < -0.29) ? 101.815 :
                (x < -1e-11) ? 127.0 :
                (x < -1e-96) ? 181.0 :
                317.0;

    double zl = LogFast(-x);
    double t = -1.0 - zl;
    double ts = sqrt(t);
    double approx = zl - (2.0 * ts) / (s2 + (c13 - t
        / (270.0 + ts * a)) * ts);

    return approx;
#else
    // Compute t
    double t = -1.0 - LogFast(-x);

    // Compute x1
    static constexpr double P1[] = {
        3.00735538504012242e+00,
        -1.11079558616957819e-02,
        1.64300407704260418e-05,
        -2.87179186820882765e-08,
        2.93030279922565400e-11,
        -1.23057197181254966e-14
    };

    static constexpr double P2[] = {
        3.00154719163439900e+00,
        -8.18302417096426201e-03,
        -1.59272314739640301e-04
    };

    double x1;
    if (t < 13.0)
    {
        x1 = P2[2];
        for (size_t i = 0; i < 2; i++)
            x1 = x1 * t + P2[1 - i];
    }
    else
    {
        x1 = P1[5];
        for (size_t i = 0; i < 5; i++)
            x1 = x1 * t + P1[4 - i];
    }

    x1 = 1.0 / x1;

    double approx = sqrt(t * 0.5);
    approx = approx * x1 + 1.0;
    approx = 6.0 / approx;
    approx = approx - t;
    approx = approx - 7.0;

    return approx;

#endif
}

double GeneralWm1(double x)
{
    double w = Approx(x);

    // === Fritsch Iteration ===
    static constexpr double c23 = 2.0 / 3.0;
    static constexpr double smallScale = 4611686018427387904.0;     // 2^62
    static constexpr double smallOffset = 42.975125194716609184;    // ln(2^62)

    bool isSmall = x > -1e-300;
    if (isSmall)
        x *= smallScale;
    double zn = log(x / w);
    if (isSmall)
        zn -= smallOffset;
    zn -= w;

    double temp = 1.0 + w;
    double temp2 = temp + c23 * zn;
    temp2 = 2.0 * temp * temp2;
    w = w * (1.0 + (zn / temp) * (temp2 - zn) / (temp2 - 2.0 * zn));

    return w;
}

double AddEm(double x)
{
    static constexpr double emHigh = 0.36787944117144232160;
    static constexpr double emLow = -1.2428753672788363168e-17;

    return (x + emHigh) + emLow;
}

double NearBranchSeries(double p)
{
    static constexpr double P[] = {
        -1,
        -1.0000000000000002,
        -0.333333333333167,
        -0.15277777780222743,
        -0.07962962816509019,
        -0.04450236076137409,
        -0.02598385679312408,
        -0.015645891525841612,
        -0.009535109600029509,
        -0.006457155334744418,
        -0.0021889836393367036,
        -0.006358886656862373,
        0.004232489664240261,
        -0.005427367255942878
    };

    // Evaluate polynomial using Horner's Method
    double value = P[13];
    for (size_t i = 0; i < 13; i++)
        value = value * p + P[12 - i];

    return value;
}

double NearBranchWm1(double x)
{
    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    return NearBranchSeries(p);
}

double MuirWm1(double x)
{
    return (x < -0.346) ? NearBranchWm1(x) : GeneralWm1(x);
}
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
        1.90672114624406,
        8.128607188367015,
        9.353551895034446,
        2.328717384810457,
        0.02574163263511813
    };

    static constexpr double Q[] = {
        1.9068870342309603,
        10.0360068543767,
        16.522152695631256,
        8.873280807877627,
        1
    };

    double numer = P[5];
    for (size_t i = 0; i < 5; i++)
        numer = numer * x + P[4 - i];

    double denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * x + Q[3 - i];

    return numer / denom;
}

static inline double SecondApprox(double x)
{
    static constexpr double P[] = {
        266.74662101711755,
        180.72015154289477,
        81.03709502548347,
        0.9987349881680496
    };
    static constexpr double Q[] = {
        438.5489337661638,
        87.11265453382124,
        1
    };

    double logX = log(x);

    double numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * logX + P[2 - i];

    double denom = Q[2];
    for (size_t i = 0; i < 2; i++)
        denom = denom * logX + Q[1 - i];

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
        0.999999999999971,
        -0.3333333333290947,
        0.1527777775369985,
        -0.07962962252727342,
        0.04450218977764427,
        -0.025983294377123514,
        0.015624730774750031,
        -0.009558915212389768,
        0.005798644499473038,
        -0.0032495150699457884,
        0.0014324144431266657,
        -0.00034866137251673187
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

    return (x < -0.342) ? NearBranchW0(x) : GeneralW0(x);
}
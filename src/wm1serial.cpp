#include <cmath>

static inline double Approx(double x)
{
    // === Constants ===
    static constexpr double P[] = {
        -3.8275093535321854,
        -6.416670530095372,
        -3.9438756505405625,
        -0.9985968226777492
    };
    static constexpr double Q = 3.827222653404185;
    // =================

    double t = sqrt(-1 - log(-x));
    double w = P[3];
    for (size_t i = 0; i < 3; i++)
        w = w * t + P[2 - i];

    return w / (t + Q);
}

static inline double GeneralWm1(double x)
{
    double w = Approx(x);

    // === Fritsch Iteration ===
    static constexpr double c23 = 2.0 / 3.0;
    static constexpr double smallScale = 4611686018427387904.0;     // 2^62
    static constexpr double smallOffset = 42.975125194716609184;    // ln(2^62)

    double zn = (x > -1e-300) ? log(x * smallScale / w) - smallOffset - w : log(x / w) - w;
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

static inline double NearBranchWm1(double x)
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

    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    // Evaluate polynomial using Horner's Method
    double value = P[13];
    for (size_t i = 0; i < 13; i++)
        value = value * p + P[12 - i];

    return value;
}

double MuirWm1(double x)
{
    return (x < -0.346) ? NearBranchWm1(x) : GeneralWm1(x);
}
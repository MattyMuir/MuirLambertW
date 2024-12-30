#include <cmath>

static inline double AddEm(double x)
{
    static constexpr double emHigh = 0.36787944117144232160;
    static constexpr double emLow = -1.2428753672788363168e-17;

    return (x + emHigh) + emLow;
}

// (-0.3578794411714423215955237701, -0.3178794411714423215955237]
static inline double Approx1(double x)
{
    static constexpr double P[] = {
        -0.9999999966370994,
        0.9999998919599464,
        -0.33333176121261077,
        0.15276407720891572,
        -0.07954996871679511,
        0.04417646254148405,
        -0.0250210964451865,
        0.013549251909039679,
        -0.006295120386207784,
        0.0021121327899401854,
        -0.0003654774623776282
    };

    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    double value = P[10];
    for (size_t i = 0; i < 10; i++)
        value = value * p + P[9 - i];

    return value;
}

// (-0.3178794411714423215955237, -0.2]
static inline double Approx2(double x)
{
    static constexpr double P[] = {
        -0.9999996784581062,
        0.9999943491115374,
        -0.3332876382845283,
        0.1525523236702569,
        -0.07887051413761093,
        0.04265265610752642,
        -0.02260517433837137,
        0.010887793759798498,
        -0.004370059467087887,
        0.0013206508704697916,
        -0.00026039124547708805,
        2.4729425099979306e-05
    };

    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    double value = P[11];
    for (size_t i = 0; i < 11; i++)
        value = value * p + P[10 - i];

    return value;
}

// (-0.2, -0.15]
static inline double Approx3(double x)
{
    static constexpr double P[] = {
        -2.00000016065932383e-01L,
        -2.43056034919570645e+00L,
        -9.29445973019267377e+00L,
        -1.39873100223530626e+01L,
        -7.13659960146842413e+00L,
        -1.70369315166269611e-01L
    };

    static constexpr double Q[] = {
        1.00000000000000000e+00L,
        7.15279775097383673e+00L,
        1.82082236178655261e+01L,
        1.92070309700806545e+01L,
        6.87274343200473953e+00L
    };

    double numer = P[5];
    for (size_t i = 0; i < 5; i++)
        numer = numer * x + P[4 - i];

    double denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * x + Q[3 - i];

    return x * (numer / denom + 1.2);
}

// (-0.15, -0.1]
static inline double Approx4(double x)
{
    static constexpr double P[] = {
        -1.00000000204783974e-01L,
        -1.68062750181694521e+00L,
        -6.93814797096171538e+00L,
        -1.03738544312384574e+01L,
        -4.96135934908649410e+00L,
        -8.96693172288665332e-02L
    };

    static constexpr double Q[] = {
        1.00000000000000000e+00L,
        6.80627486726167287e+00L,
        1.63187261601690365e+01L,
        1.59786488758684961e+01L,
        5.18962689407092801e+00L
    };

    double numer = P[5];
    for (size_t i = 0; i < 5; i++)
        numer = numer * x + P[4 - i];

    double denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * x + Q[3 - i];

    return x * (numer / denom + 1.1);
}

// (-0.1, -0.051)
static inline double Approx5(double x)
{
    static constexpr double Y = 1.08633995056152344e+00;
    static constexpr double P[] = {
        -8.63399505615014331e-02,
        -1.64303871814816464e+00,
        -7.71247913918273738e+00,
        -1.41014495545382454e+01,
        -1.02269079949257616e+01,
        -2.17236002836306691e+00,
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        7.44775406945739243e+00,
        2.04392643087266541e+01,
        2.51001961077774193e+01,
        1.31256080849023319e+01,
        2.11640324843601588e+00,
    };

    double numer = P[5];
    for (size_t i = 0; i < 5; i++)
        numer = numer * x + P[4 - i];

    double denom = Q[5];
    for (size_t i = 0; i < 5; i++)
        denom = denom * x + Q[4 - i];

    return x * (Y + numer / denom);
}

// [-0.051, 0.051]
static inline double Approx6(double x)
{
	double result =
        x * (1. - // j1 x^1
        x * (1. -  // j2 x^2
        x * (1.5 - // 3/2 // j3 x^3
        x * (2.6666666666666666667 -  // 8/3 // j4
        x * (5.2083333333333333333 - // -125/24 // j5
        x * (10.8 - // j6
        x * (23.343055555555555556 - // j7
        x * (52.012698412698412698 - // j8
        x * (118.62522321428571429 - // j9
        x * (275.57319223985890653 - // j10
        x * (649.78717234347442681 - // j11
        x * (1551.1605194805194805 - // j12
        x * (3741.4497029592385495 - // j13
        x * (9104.5002411580189358 - // j14
        x * (22324.308512706601434 - // j15
        x * (55103.621972903835338 - // j16
        x * 136808.86090394293563)))))))))))))))); // j17 x^17

    return result;
}

// (0.051, 0.5)
static inline double Approx7(double x)
{
    static constexpr double Y = 8.19659233093261719e-01;
    static constexpr double P[] = {
        1.80340766906685177e-01,
        3.28178241493119307e-01,
        -2.19153620687139706e+00,
        -7.24750929074563990e+00,
        -7.28395876262524204e+00,
        -2.57417169492512916e+00,
        -2.31606948888704503e-01
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        7.36482529307436604e+00,
        2.03686007856430677e+01,
        2.62864592096657307e+01,
        1.59742041380858333e+01,
        4.03760534788374589e+00,
        2.91327346750475362e-01
    };

    double numer = P[6];
    for (size_t i = 0; i < 6; i++)
        numer = numer * x + P[5 - i];

    double denom = Q[6];
    for (size_t i = 0; i < 6; i++)
        denom = denom * x + Q[5 - i];

    return x * (Y + numer / denom);
}

// [0.5, 2)
static inline double Approx8(double x)
{
    static constexpr double Y = 5.50335884094238281e-01;
    static constexpr double P[] = {
        4.49664083944098322e-01,
        1.90417666196776909e+00,
        1.99951368798255994e+00,
        -6.91217310299270265e-01,
        -1.88533935998617058e+00,
        -7.96743968047750836e-01,
        -1.02891726031055254e-01,
        -3.09156013592636568e-03
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        6.45854489419584014e+00,
        1.54739232422116048e+01,
        1.72606164253337843e+01,
        9.29427055609544096e+00,
        2.29040824649748117e+00,
        2.21610620995418981e-01,
        5.70597669908194213e-03
    };

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * x + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * x + Q[6 - i];

    return x * (Y + numer / denom);
}

// [2, 6)
static inline double Approx9(double x)
{
    static constexpr double Y = 1.16239356994628906e+00;
    static constexpr double P[] = {
        -1.16230494982099475e+00,
        -3.38528144432561136e+00,
        -2.55653717293161565e+00,
        -3.06755172989214189e-01,
        1.73149743765268289e-01,
        3.76906042860014206e-02,
        1.84552217624706666e-03,
        1.69434126904822116e-05,
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        3.77187616711220819e+00,
        4.58799960260143701e+00,
        2.24101228462292447e+00,
        4.54794195426212385e-01,
        3.60761772095963982e-02,
        9.25176499518388571e-04,
        4.43611344705509378e-06,
    };

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * x + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * x + Q[6 - i];

    return Y + numer / denom;
}

// [6, 18)
static inline double Approx10(double x)
{
    static constexpr double Y = 1.80937194824218750e+00;
    static constexpr double P[] = {
        -1.80690935424793635e+00,
        -3.66995929380314602e+00,
        -1.93842957940149781e+00,
        -2.94269984375794040e-01,
        1.81224710627677778e-03,
        2.48166798603547447e-03,
        1.15806592415397245e-04,
        1.43105573216815533e-06,
        3.47281483428369604e-09
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        2.57319080723908597e+00,
        1.96724528442680658e+00,
        5.84501352882650722e-01,
        7.37152837939206240e-02,
        3.97368430940416778e-03,
        8.54941838187085088e-05,
        6.05713225608426678e-07,
        8.17517283816615732e-10
    };

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * x + P[7 - i];

    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
        denom = denom * x + Q[7 - i];

    return Y + numer / denom;
}

// [18, 9897.12905874)
static inline double Approx11(double x)
{
    static constexpr double Y = -1.40297317504882812e+00;
    static constexpr double P[] = {
        1.97011826279311924e+00,
        1.05639945701546704e+00,
        3.33434529073196304e-01,
        3.34619153200386816e-02,
        -5.36238353781326675e-03,
        -2.43901294871308604e-03,
        -2.13762095619085404e-04,
        -4.85531936495542274e-06,
        -2.02473518491905386e-08,
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        8.60107275833921618e-01,
        4.10420467985504373e-01,
        1.18444884081994841e-01,
        2.16966505556021046e-02,
        2.24529766630769097e-03,
        9.82045090226437614e-05,
        1.36363515125489502e-06,
        3.44200749053237945e-09,
    };

    double t = log(x);

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * t + P[7 - i];

    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
        denom = denom * t + Q[7 - i];

    return t + Y + numer / denom;
}

// [9897.12905874, 7.896296e+13)
static inline double Approx12(double x)
{
    static constexpr double Y = -2.73572921752929688e+00;
    static constexpr double P[] = {
        3.30547638424076217e+00,
        1.64050071277550167e+00,
        4.57149576470736039e-01,
        4.03821227745424840e-02,
        -4.99664976882514362e-04,
        -1.28527893803052956e-04,
        -2.95470325373338738e-06,
        -1.76662025550202762e-08,
        -1.98721972463709290e-11,
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        6.91472559412458759e-01,
        2.48154578891676774e-01,
        4.60893578284335263e-02,
        3.60207838982301946e-03,
        1.13001153242430471e-04,
        1.33690948263488455e-06,
        4.97253225968548872e-09,
        3.39460723731970550e-12,
    };

    double t = log(x);

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * t + P[7 - i];

    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
        denom = denom * t + Q[7 - i];

    return t + Y + numer / denom;
}

// [7.896296e+13, 2.6881171e+43)
static inline double Approx13(double x)
{
    static constexpr double Y = -4.01286315917968750e+00;
    static constexpr double P[] = {
        5.07714858354309672e+00,
        -3.32994414518701458e+00,
        -8.61170416909864451e-01,
        -4.01139705309486142e-02,
        -1.85374201771834585e-04,
        1.08824145844270666e-05,
        1.17216905810452396e-07,
        2.97998248101385990e-10,
        1.42294856434176682e-13,
    };
    static constexpr double Q[] = {
        1.00000000000000000e+00,
        -4.85840770639861485e-01,
        -3.18714850604827580e-01,
        -3.20966129264610534e-02,
        -1.06276178044267895e-03,
        -1.33597828642644955e-05,
        -6.27900905346219472e-08,
        -9.35271498075378319e-11,
        -2.60648331090076845e-14,
    };

    double t = log(x);

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * t + P[7 - i];

    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
        denom = denom * t + Q[7 - i];

    return t + Y + numer / denom;
}

// [2.6881171e+43, DBL_MAX)
static inline double Approx14(double x)
{
    static constexpr double Y = -5.70115661621093750e+00;
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
        -4.06109208815303157e-20,
        -2.21552699006496737e-24,
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
        1.11775518708172009e-20,
        3.78250395617836059e-25,
    };

    double t = log(x);

    double numer = P[10];
    for (size_t i = 0; i < 10; i++)
        numer = numer * t + P[9 - i];

    double denom = Q[10];
    for (size_t i = 0; i < 10; i++)
        denom = denom * t + Q[9 - i];

    return t + Y + numer / denom;
}

static inline double NearBranchW0(double x)
{
    static constexpr double P[] = {
        -0.9999999999999999,
        0.9999999999999364,
        -0.33333333332245435,
        0.15277777704999906,
        -0.07962960455991744,
        0.044501808551933024,
        -0.02597829873566169,
        0.015582883821839718,
        -0.00933353422702553,
        0.00503749851616291,
        -0.0017510690299764934
    };

    static constexpr double s2e = 2.331643981597124;
    double p = sqrt(AddEm(x)) * s2e;

    double value = P[10];
    for (size_t i = 0; i < 10; i++)
        value = value * p + P[9 - i];

    return value;
}

double MuirW0(double x)
{
    if (x < -0.3578794411714423215955237701)
        return NearBranchW0(x);

    if (x < 0.5)
    {
        if (x < -0.1)
        {
            if (x < -0.2)
            {
                if (x < -0.3178794411714423215955237) return Approx1(x);
                return Approx2(x);
            }
            if (x < -0.15) return Approx3(x);
            return Approx4(x);
        }
        if (x < 0.051)
        {
            if (x < -0.051) return Approx5(x);
            return Approx6(x);
        }
        return Approx7(x);
    }
    if (x < 9897.12905874)
    {
        if (x < 6)
        {
            if (x < 2) return Approx8(x);
            return Approx9(x);
        }
        if (x < 18) return Approx10(x);
        return Approx11(x);
    }
    if (x < 2.6881171e+43)
    {
        if (x < 7.896296e+13) return Approx12(x);
        return Approx13(x);
    }
    return Approx14(x);
}
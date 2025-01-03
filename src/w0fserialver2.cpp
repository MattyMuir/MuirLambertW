#include <cmath>

static inline float AddEm(float x)
{
    static constexpr float emHigh = 0.36787945f;
    static constexpr float emLow = -9.149756e-09f;

    return (x + emHigh) + emLow;
}

// (-0.3578794411714423215955237701, -0.27]
static inline float Approx1(float x)
{
    static constexpr float P[] = {
        -0.9999900179300805f,
        0.999825992034222f,
        -0.3320510417108332f,
        0.1475469965715199f,
        -0.06661836087728089f,
        0.023903508735882057f,
        -0.004588270097285815f
    };

    static constexpr float s2e = 2.331644f;
    float p = sqrtf(AddEm(x)) * s2e;

    float value = P[6];
    for (size_t i = 0; i < 6; i++)
        value = value * p + P[5 - i];

    return value;
}

// (-0.27, -0.051)
static inline float Approx2(float x)
{
    static constexpr float Y = 1.255809784e+00f;
    static constexpr float P[] = {
        -2.558083412e-01f,
        -2.306524098e+00f,
        -5.630887033e+00f,
        -3.803974556e+00f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        5.107680783e+00f,
        7.914062868e+00f,
        3.501498501e+00f,
    };

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * x + P[2 - i];

    float denom = Q[3];
    for (size_t i = 0; i < 3; i++)
        denom = denom * x + Q[2 - i];

    return x * (Y + numer / denom);
}

// [-0.051, 0.05]
static inline float Approx3(float x)
{
    static constexpr float P[] = {
        0.0f,
        1.0f,
        -1.0000000424879305f,
        1.499951592926675f,
        -2.666387307065152f,
        5.283832465824341f,
        -11.050509163760704f
    };

    float numer = P[6];
    for (size_t i = 0; i < 6; i++)
        numer = numer * x + P[5 - i];

    return numer;
}

// (0.05, 0.5)
static inline float Approx4(float x)
{
    static constexpr float Y = 8.196592331e-01f;
    static constexpr float P[] = {
        1.803388345e-01f,
        -4.820256838e-01f,
        -1.068349741e+00f,
        -3.506624319e-02f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        2.871703469e+00f,
        1.690949264e+00f,
    };

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * x + P[2 - i];

    float denom = Q[2];
    for (size_t i = 0; i < 2; i++)
        denom = denom * x + Q[1 - i];

    return x * (Y + numer / denom);
}

// [0.5, 2)
static inline float Approx5(float x)
{
    static constexpr float Y = 5.503368378e-01f;
    static constexpr float P[] = {
        4.493332766e-01f,
        2.543432707e-01f,
        -4.808788799e-01f,
        -1.244425316e-01f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        2.780661241e+00f,
        1.830840318e+00f,
        2.407221031e-01f,
    };

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * x + P[2 - i];

    float denom = Q[3];
    for (size_t i = 0; i < 3; i++)
        denom = denom * x + Q[2 - i];

    return x * (Y + numer / denom);
}

// [2, 6)
static inline float Approx6(float x)
{
    static constexpr float Y = 1.162393570e+00f;
    static constexpr float P[] = {
        -1.144183394e+00f,
        -4.712732855e-01f,
        1.563162512e-01f,
        1.434010911e-02f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        1.192626340e+00f,
        2.295580708e-01f,
        5.477869455e-03f,
    };

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * x + P[2 - i];

    float denom = Q[3];
    for (size_t i = 0; i < 3; i++)
        denom = denom * x + Q[2 - i];

    return Y + numer / denom;
}

// [6, 18)
static inline float Approx7(float x)
{
    static constexpr float Y = 1.809371948e+00f;
    static constexpr float P[] = {
        -1.689291769e+00f,
        -3.337812742e-01f,
        3.151434873e-02f,
        1.134178734e-03f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        5.716915685e-01f,
        4.489521292e-02f,
        4.076716763e-04f,
    };

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * x + P[2 - i];

    float denom = Q[3];
    for (size_t i = 0; i < 3; i++)
        denom = denom * x + Q[2 - i];

    return Y + numer / denom;
}

// [18, 9897.12905874)
static inline float Approx8(float x)
{
    static constexpr float Y = -1.402973175e+00f;
    static constexpr float P[] = {
        1.966174312e+00f,
        2.350864728e-01f,
        -5.098074353e-02f,
        -1.054818339e-02f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        4.388208264e-01f,
        8.316639634e-02f,
        3.397187918e-03f,
        -1.321489743e-05f,
    };

    float t = logf(x);

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * t + P[2 - i];

    float denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * t + Q[3 - i];

    return t + Y + numer / denom;
}

// [9897.12905874, 7.896296e+13)
static inline float Approx9(float x)
{
    static constexpr float Y = -2.735729218e+00f;
    static constexpr float P[] = {
        3.424903470e+00f,
        7.525631787e-02f,
        -1.427309584e-02f,
        -1.435974178e-05f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        2.514005579e-01f,
        6.118994652e-03f,
        -1.357889535e-05f,
        7.312865624e-08f,
    };

    float t = logf(x);

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * t + P[2 - i];

    float denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * t + Q[3 - i];

    return t + Y + numer / denom;
}

// [7.896296e+13, FLT_MAX)
static inline float Approx10(float x)
{
    static constexpr float Y = -4.012863159e+00f;
    static constexpr float P[] = {
        4.431629226e+00f,
        2.756690487e-01f,
        -2.992956930e-03f,
        -4.912259384e-05f,
    };
    static constexpr float Q[] = {
        1.000000000e+00f,
        2.015434591e-01f,
        4.949426142e-03f,
        1.609659944e-05f,
        -5.111523436e-09f,
    };

    float t = logf(x);

    float numer = P[3];
    for (size_t i = 0; i < 3; i++)
        numer = numer * t + P[2 - i];

    float denom = Q[4];
    for (size_t i = 0; i < 4; i++)
        denom = denom * t + Q[3 - i];

    return t + Y + numer / denom;
}

static inline float NearBranchW0(float x)
{
    static constexpr float P[] = {
        -0.9999999541597892f,
        0.9999903356903502f,
        -0.33300152298966174f,
        0.14870683416255348f,
        -0.0584015853073767f
    };

    static constexpr float s2e = 2.331644f;
    float p = sqrtf(AddEm(x)) * s2e;

    float value = P[4];
    for (size_t i = 0; i < 4; i++)
        value = value * p + P[3 - i];

    return value;
}

float MuirW0v2(float x)
{
    if (x < -0.3578794411714423215955237701f)
        return NearBranchW0(x);

    if (x < 2.0f)
    {
        if (x < 0.05f)
        {
            if (x < -0.051f)
            {
                if (x < -0.27f) return Approx1(x);
                return Approx2(x);
            }
            return Approx3(x);
        }
        if (x < 0.5f) return Approx4(x);
        return Approx5(x);
    }
    if (x < 9897.12905874f)
    {
        if (x < 18.0f)
        {
            if (x < 6.0f) return Approx6(x);
            return Approx7(x);
        }
        return Approx8(x);
    }
    if (x < 7.896296e+13f) return Approx9(x);
    return Approx10(x);
}
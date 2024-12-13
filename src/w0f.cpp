#include <immintrin.h>

struct double8
{
    __m256d l, h;
};

static inline double8 Set1(double v)
{
    return { _mm256_set1_pd(v), _mm256_set1_pd(v) };
}

static inline double8 fmadd(double8 x, double8 y, double8 z)
{
    return { _mm256_fmadd_pd(x.l, y.l, z.l), _mm256_fmadd_pd(x.h, y.h, z.h) };
}

static inline double8 FloatToDouble(__m256 p)
{
    __m256d l = _mm256_cvtps_pd(_mm256_castps256_ps128(p));
    __m256d h = _mm256_cvtps_pd(_mm256_extractf128_ps(p, 1));
    return { l, h };
}

static inline __m256 DoubleToFloat(double8 p)
{
    return _mm256_set_m128(_mm256_cvtpd_ps(p.h), _mm256_cvtpd_ps(p.l));
}
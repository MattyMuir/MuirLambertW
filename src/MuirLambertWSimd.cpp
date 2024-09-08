#include "MuirLambertWSimd.h"

#include <cstdint>
#include <iostream>
#include <bit>
#include <limits>
#include <immintrin.h>

#include "mysleef.h"

#include "io.h"

#define LESS 0x11
#define EQUAL 0x0
#define GREATER 0x1E

#define ANY(p) (bool)(_mm256_movemask_pd(p))
#define ALL_ONE _mm256_cmp_pd(_mm256_setzero_pd(), _mm256_setzero_pd(), 0x0)
#define NOT(p) _mm256_andnot_pd(p, ALL_ONE)

#define SHUFFLE_INT(a, b, i) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), i))

static inline __m256i Lzcnt(__m256i a)
{
    uint64_t c0 = std::countl_zero((uint64_t)_mm256_extract_epi64(a, 0));
    uint64_t c1 = std::countl_zero((uint64_t)_mm256_extract_epi64(a, 1));
    uint64_t c2 = std::countl_zero((uint64_t)_mm256_extract_epi64(a, 2));
    uint64_t c3 = std::countl_zero((uint64_t)_mm256_extract_epi64(a, 3));

    return _mm256_setr_epi64x(c0, c1, c2, c3);
}

static inline __m256i Mantissa(__m256d x)
{
    static constexpr uint64_t MantissaMask = 0xfffffffffffff;
    static constexpr uint64_t ImpliedOne = (1ULL << 52);

    __m256i punn = _mm256_castpd_si256(x);
    punn = _mm256_and_si256(punn, _mm256_set1_epi64x(MantissaMask));
    punn = _mm256_or_si256(punn, _mm256_set1_epi64x(ImpliedOne));

    return punn;
}

static inline __m256d AddEm(__m256d x)
{
    // Constants
    static constexpr uint64_t emLow = 0xc6aeb7b1e0a4153e;
    static constexpr uint64_t emHigh = 0x00178b56362cef37;

    // Subtract mantissas
    __m256i fullLow = _mm256_set1_epi64x(emLow);
    __m256i fullHigh = _mm256_sub_epi64(_mm256_set1_epi64x(emHigh), Mantissa(x));

    // Get shift amount
    __m256i shamt = _mm256_sub_epi64(Lzcnt(fullHigh), _mm256_set1_epi64x(11));

    // Perform shift
    __m256i newMan = _mm256_sllv_epi64(fullHigh, shamt);
    __m256i oppShamt = _mm256_sub_epi64(_mm256_set1_epi64x(64), shamt);
    newMan = _mm256_or_si256(newMan, _mm256_srlv_epi64(fullLow, oppShamt));

    // New exponent
    __m256i biasedExp = _mm256_sub_epi64(_mm256_set1_epi64x(1021), shamt);

    // Combine mantissa and exponent
    __m256i result = _mm256_and_si256(newMan, _mm256_set1_epi64x(0xfffffffffffff));
    result = _mm256_or_si256(result, _mm256_slli_epi64(biasedExp, 52));

    return _mm256_castsi256_pd(result);
}

template <size_t MaxOrder>
static __m256d NearBranchSeries(__m256d p)
{
    // === Constants ===
    static constexpr double q[] = {
      -1,
      +1,
      -0.333333333333333333,
      +0.152777777777777778,
      -0.0796296296296296296,
      +0.0445023148148148148,
      -0.0259847148736037625,
      +0.0156356325323339212,
      -0.00961689202429943171,
      +0.00601454325295611786,
      -0.00381129803489199923,
      +0.00244087799114398267,
      -0.00157693034468678425,
      +0.00102626332050760715,
      -0.000672061631156136204,
      +0.000442473061814620910,
      -0.000292677224729627445,
      +0.000194387276054539318,
      -0.000129574266852748819,
      +0.0000866503580520812717,
      -0.0000581136075044138168,
      +0.00003907668486743905163539558,
      -0.00002633806474723109873858408,
      +0.00001779034580507958540073628,
      -0.00001204035273955997694227412
    };
    // =================

    // Evaluate polynomial using Horner's Method
    __m256d value = _mm256_fmadd_pd(p, _mm256_set1_pd(q[MaxOrder]), _mm256_set1_pd(q[MaxOrder - 1]));

    for (int64_t i = MaxOrder - 2; i >= 0; i--)
        value = _mm256_fmadd_pd(p, value, _mm256_set1_pd(q[i]));

    return value;
}

static __m256d NearBranchW0(__m256d x)
{
    // Compute p = sqrt(2 e x + 2)
    static constexpr double s2e = 2.331643981597124;
    __m256d p = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(s2e));

    return NearBranchSeries<15>(p);
}

static inline __m256d Epi64ToPd(__m256i x)
{
    __m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    x = _mm256_permutevar8x32_epi32(x, perm);
    return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

// Calculates an approximation to ln(x) for normalized doubles
static __m256d LogApprox(__m256d x)
{
    __m256i punn = _mm256_castpd_si256(x);
    __m128i punnTop = _mm256_extractf128_si256(punn, 1);
    __m128i compressed = _mm256_castsi256_si128(SHUFFLE_INT(punn, _mm256_castsi128_si256(punnTop), 221));
    compressed = _mm_sub_epi32(compressed, _mm_set1_epi32(1072632447));
    return _mm256_mul_pd(_mm256_cvtepi32_pd(compressed), _mm256_set1_pd(6.610368362777016e-7));
}

static __m256d BetterLogApprox(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    __m256d c = _mm256_set1_pd(-0.238626940484713889774143051645); // This constant ensures the approximation is exact for ln(1/e)
    __m256d one = _mm256_set1_pd(1.0);
    __m256d two = _mm256_set1_pd(2.0);
    // =================

    // Extract exponent
    __m256i punn = _mm256_castpd_si256(x);
    __m256d exp = Epi64ToPd(_mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023)));

    // Extract mantissa
    punn = _mm256_and_si256(punn, _mm256_set1_epi64x(0x3FFFFFFFFFFFFFFF));
    punn = _mm256_or_si256(punn, _mm256_set1_epi64x(0x3FF0000000000000));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Compute approximation
    __m256d res = _mm256_sub_pd(_mm256_add_pd(exp, mantissa), _mm256_set1_pd(1.0));
    res = _mm256_mul_pd(res, ln2);

    __m256d mSubOne = _mm256_sub_pd(mantissa, one);
    __m256d temp = _mm256_mul_pd(mSubOne, _mm256_sub_pd(mantissa, two));
    res = _mm256_fmadd_pd(temp, c, res);
    return res;
}

__m256d W0Iterations(__m256d x, __m256d w)
{
    // Constants
    static constexpr double c23 = 0.6666666666666666;
    __m256d one = _mm256_set1_pd(1.0);

    // === Fritsch Iteration ===
    __m256d xov = _mm256_div_pd(x, w);
    __m256d valsEq = _mm256_cmp_pd(x, w, EQUAL);
    xov = _mm256_blendv_pd(xov, one, valsEq);
    __m256d zn = _mm256_sub_pd(Sleef_logd4_u35avx2(xov), w);

    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, _mm256_set1_pd(c23), temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    // === Halley Iteration ===
    __m256d expW = Sleef_expd4_u10avx2(w);
    __m256d wExpWSubX = _mm256_fmsub_pd(w, expW, x);

    __m256d newW = _mm256_add_pd(w, _mm256_set1_pd(2));
    newW = _mm256_mul_pd(newW, wExpWSubX);
    newW = _mm256_div_pd(newW, _mm256_fmadd_pd(w, _mm256_set1_pd(2), _mm256_set1_pd(2)));
    newW = _mm256_sub_pd(_mm256_fmadd_pd(w, expW, expW), newW);
    newW = _mm256_div_pd(wExpWSubX, newW);
    newW = _mm256_sub_pd(w, newW);
    w = newW;

    return w;
}

__m256d WM1Iterations(__m256d x, __m256d w)
{
    // Constants
    static constexpr double c23 = 0.6666666666666666;
    __m256d one = _mm256_set1_pd(1.0);

    // === Fritsch Iteration ===
    __m256d xov = _mm256_div_pd(x, w);
    __m256d valsEq = _mm256_cmp_pd(x, w, EQUAL);
    xov = _mm256_blendv_pd(xov, one, valsEq);
    __m256d zn = _mm256_sub_pd(Sleef_logd4_u35avx2(xov), w); // For denormalized arguments, first multiply by a fixed factor and subtract after
    
    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, _mm256_set1_pd(c23), temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    return w;
}

static __m256d GeneralW0(__m256d x)
{
    // === Constants ===
    static constexpr double e = 2.7182818284590452353602874;
    static constexpr double em = 0.36787944117144232160;
    static constexpr double rt2e = 2.3316439815971242034;
    static constexpr double a = -2.1;
    static constexpr double b = -0.21846942280995411989;
    static constexpr double c = -0.92211551407955841243;
    static constexpr double d = -1.4142135623730950488;
    __m256d one = _mm256_set1_pd(1.0);
    // =================

    // First approximation
    __m256d lx = LogApprox(_mm256_add_pd(x, _mm256_set1_pd(e)));
    __m256d scale = _mm256_blendv_pd(_mm256_set1_pd(-1), _mm256_set1_pd(-0.63), _mm256_cmp_pd(x, _mm256_set1_pd(100), LESS));
    __m256d approx1 = _mm256_fmadd_pd(LogApprox(lx), scale, lx);
    approx1 = _mm256_mul_pd(_mm256_div_pd(x, _mm256_add_pd(x, one)), approx1);

    // Second approximation
    __m256d xpem = _mm256_add_pd(x, _mm256_set1_pd(em));
    __m256d reta = _mm256_mul_pd(_mm256_sqrt_pd(xpem), _mm256_set1_pd(rt2e));
    __m256d approx2 = _mm256_fmadd_pd(x, _mm256_set1_pd(a), _mm256_set1_pd(b));
    approx2 = _mm256_fmadd_pd(x, approx2, _mm256_set1_pd(c));
    approx2 = _mm256_fmadd_pd(x, approx2, _mm256_set1_pd(d));
    approx2 = _mm256_add_pd(approx2, reta);

    // Blend approximations
    __m256d useLarge = _mm256_cmp_pd(x, _mm256_set1_pd(-0.185), GREATER);
    __m256d value = _mm256_blendv_pd(approx2, approx1, useLarge);

    // Iteration
    return W0Iterations(x, value);
}

__m256d MuirLambertW0Simd(__m256d x)
{
    // === Constants ===
    const __m256d infPack = _mm256_set1_pd(std::numeric_limits<double>::infinity());

    // Get masks for values near branch
    __m256d notInf = _mm256_cmp_pd(x, infPack, LESS);
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.362), LESS);
    isNearBranch = _mm256_and_pd(isNearBranch, notInf);
    __m256d isNotNearBranch = _mm256_andnot_pd(isNearBranch, notInf);

    // Compute values
    __m256d values, nbValues;
    if (ANY(isNearBranch))
        nbValues = NearBranchW0(x);
    if (ANY(isNotNearBranch))
        values = GeneralW0(x);

    // Blend
    __m256d results = _mm256_blendv_pd(values, nbValues, isNearBranch);
    results = _mm256_blendv_pd(infPack, results, notInf);

    return results;
}

static __m256d NearBranchWM1(__m256d x)
{
    static constexpr double s2e = 2.331643981597124;

    // Compute p = -sqrt(2 e x + 2)
    __m256d p = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(-s2e));

    return NearBranchSeries<22>(p);
}

__m256d GeneralWM1(__m256d x)
{
    // Constants
    __m256d one = _mm256_set1_pd(1.0);
    __m256d negOne = _mm256_set1_pd(-1.0);
    __m256d c = _mm256_set1_pd(0.25);
    __m256d tonc = _mm256_set1_pd(-8);

    // Approximation
    // Determine value for 'a'
    __m256d xQuiteSmall = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-11), GREATER);
    __m256d xVerySmall = _mm256_cmp_pd(x, _mm256_set1_pd(-1e-96), GREATER);
    __m256d a = _mm256_blendv_pd(_mm256_set1_pd(127), _mm256_set1_pd(181), xQuiteSmall);
    a = _mm256_blendv_pd(a, _mm256_set1_pd(317), xVerySmall);

    __m256d logX = Sleef_logd4_u35avx2(_mm256_sub_pd(_mm256_setzero_pd(), x));
    __m256d t = _mm256_sub_pd(negOne, logX);
    __m256d sqrtT = _mm256_sqrt_pd(t);
    __m256d approx = _mm256_fmadd_pd(a, sqrtT, _mm256_set1_pd(270));
    approx = _mm256_div_pd(t, approx);
    approx = _mm256_sub_pd(_mm256_set1_pd(1.0 / 3.0), approx);
    approx = _mm256_fmadd_pd(approx, sqrtT, _mm256_set1_pd(1.4142135623730950488));
    approx = _mm256_div_pd(sqrtT, approx);
    approx = _mm256_add_pd(approx, approx);
    approx = _mm256_sub_pd(logX, approx);

    return WM1Iterations(x, approx);
}

__m256d MuirLambertWM1Simd(__m256d x)
{
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.345), LESS);
    
    // Compute values
    __m256d values, nbValues;
    if (ANY(isNearBranch))
        nbValues = NearBranchWM1(x);
    if (ANY(NOT(isNearBranch)))
        values = GeneralWM1(x);

    // Blend
    __m256d results = _mm256_blendv_pd(values, nbValues, isNearBranch);

    return results;
}
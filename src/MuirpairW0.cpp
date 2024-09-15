#include "MuirpairW0.h"

#include <cstdint>
#include <limits>

#define EQUAL 0x0
#define LESS 0x11
#define GREATER 0x1E

#define SHUFFLE_INT(a, b, i) _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), i))

// ========== Log Functions ==========
static inline __m256d Epi64ToPd(__m256i x)
{
    __m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    x = _mm256_permutevar8x32_epi32(x, perm);
    return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

static __m256d LogFast(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    // =================

    // Extract exponent
    __m256i punn = _mm256_castpd_si256(x);
    __m256d exp = Epi64ToPd(_mm256_sub_epi64(_mm256_srli_epi64(punn, 52), _mm256_set1_epi64x(1023)));

    // Extract mantissa
    punn = _mm256_and_si256(punn, _mm256_set1_epi64x(0x3FFFFFFFFFFFFFFF));
    punn = _mm256_or_si256(punn, _mm256_set1_epi64x(0x3FF0000000000000));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Compute approximation
    static constexpr double P[] = {
        -2.001981174757893,
        3.7491905722285463,
        -2.77383516599771,
        1.3494246360270206,
        -0.36414336579589207,
        0.04134449829706258
    };

    __m256d approx = _mm256_set1_pd(P[5]);
    for (size_t i = 0; i < 5; i++)
        approx = _mm256_fmadd_pd(approx, mantissa, _mm256_set1_pd(P[4 - i]));

    return _mm256_fmadd_pd(exp, ln2, approx);
}

static __m256d LogAccurate(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    __m256d one = _mm256_set1_pd(1.0);
    // =================

    // Extract exponent
    __m256d xd = _mm256_mul_pd(x, _mm256_set1_pd(1.0 / 0.75));
    __m256i dpunn = _mm256_castpd_si256(xd);
    __m256i exp = _mm256_sub_epi64(_mm256_srli_epi64(dpunn, 52), _mm256_set1_epi64x(1023));

    // Extract mantissa
    __m256i punn = _mm256_castpd_si256(x);
    punn = _mm256_sub_epi64(punn, _mm256_slli_epi64(exp, 52));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    __m256d t = _mm256_div_pd(_mm256_sub_pd(mantissa, one), _mm256_add_pd(mantissa, one));
    __m256d arg = _mm256_mul_pd(t, t);
    __m256d t3 = _mm256_mul_pd(t, arg);

    // Compute approximation
    static constexpr double P[] = {
        0.6666666666667778740063,
        0.399999999950799600689777,
        0.285714294746548025383248,
        0.222221366518767365905163,
        0.181863266251982985677316,
        0.152519917006351951593857,
        0.153487338491425068243146
    };

    __m256d approx = _mm256_set1_pd(P[6]);
    for (size_t i = 0; i < 6; i++)
        approx = _mm256_fmadd_pd(approx, arg, _mm256_set1_pd(P[5 - i]));
    approx = _mm256_fmadd_pd(t3, approx, _mm256_add_pd(t, t));

    return _mm256_fmadd_pd(Epi64ToPd(exp), ln2, approx);
}

// ========== General ==========
__m256d Abs(__m256d x)
{
    __m256d signMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7fff'ffff'ffff'ffff));
    return _mm256_and_pd(x, signMask);
}

__m256d FirstApprox(__m256d x)
{
    // === Constants ===
    __m256d e2 = _mm256_set1_pd(5.4365636569180904707);			// e * 2
    __m256d two = _mm256_set1_pd(2.0);							// 2
    // =================

#if 1
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

    __m256d reta = _mm256_sqrt_pd(_mm256_fmadd_pd(x, e2, two));
    //__m256d reta = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(2.3316439815971242034));

    __m256d approx = _mm256_set1_pd(P[10]);
    for (size_t i = 0; i < 10; i++)
        approx = _mm256_fmadd_pd(approx, reta, _mm256_set1_pd(P[9 - i]));

    // Use approx = x for arguments near zero
    __m256d isNearZero = _mm256_cmp_pd(Abs(x), _mm256_set1_pd(1e-4), LESS);
    approx = _mm256_blendv_pd(approx, x, isNearZero);

    return approx;
#else
    static constexpr double P[] = {
        -1.0,
        0.999941496640904044459909267686,
        -0.3318756444333983979255719987122574821115,
        0.1456506832093031611474742703649098984897,
        -0.0638525222312049667294786559068597853184,
        0.0236845051849372927044168335442009265535,
        -0.0067196890180458105495953091690353176091,
        0.0013810915789711295919434119738866684202,
        -0.0001988032504614401895903091910255966468,
        0.0000193723975900881618950055695904666209,
        -0.0000012116330584853344531112533111194551,
        0.0000000437661726042473028299211389949980,
        -0.0000000006928406708333021054312135438642
    };

    __m256d reta = _mm256_sqrt_pd(_mm256_fmadd_pd(x, e2, two));

    __m256d approx = _mm256_set1_pd(P[12]);
    for (size_t i = 0; i < 12; i++)
        approx = _mm256_fmadd_pd(approx, reta, _mm256_set1_pd(P[11 - i]));

    // Use approx = x for arguments near zero
    __m256d isNearZero = _mm256_cmp_pd(Abs(x), _mm256_set1_pd(1e-4), LESS);
    approx = _mm256_blendv_pd(approx, x, isNearZero);

    return approx;
#endif
}

__m256d SecondApprox(__m256d x)
{
    // === Constants ===
    __m256d y = _mm256_set1_pd(-5.70115661621093750e+00);
    __m256d offset = _mm256_set1_pd(-0.132);
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
       //-4.06109208815303157e-20,
       //-2.21552699006496737e-24,
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
       //1.11775518708172009e-20,
       //3.78250395617836059e-25,
    };

    __m256d logX = LogFast(_mm256_add_pd(x, offset));

    __m256d numer = _mm256_set1_pd(P[8]);
    __m256d denom = _mm256_set1_pd(Q[8]);
    for (size_t i = 0; i < 8; i++)
    {
        numer = _mm256_fmadd_pd(numer, logX, _mm256_set1_pd(P[7 - i]));
        denom = _mm256_fmadd_pd(denom, logX, _mm256_set1_pd(Q[7 - i]));
    }

    __m256d approx = _mm256_div_pd(numer, denom);
    approx = _mm256_add_pd(_mm256_add_pd(approx, y), logX);
    return approx;
}

__m256d GeneralW0(__m256d x)
{
    __m256d isOver20 = _mm256_cmp_pd(x, _mm256_set1_pd(20.0), GREATER);
    uint32_t isOver20Mask = _mm256_movemask_pd(isOver20);

    __m256d w;
    switch (isOver20Mask)
    {
    case 0b0000:
        w = FirstApprox(x);
        break;
    case 0b1111:
        w = SecondApprox(x);
        break;
    default:
        w = _mm256_blendv_pd(FirstApprox(x), SecondApprox(x), isOver20);
    }

    // Constants
    __m256d c23 = _mm256_set1_pd(2.0 / 3.0);
    __m256d one = _mm256_set1_pd(1.0);

    // === Fritsch Iteration ===
    __m256d xov = _mm256_div_pd(x, w);
    __m256d zn = _mm256_sub_pd(LogAccurate(xov), w);

    __m256d temp = _mm256_add_pd(w, one);
    __m256d temp2 = _mm256_fmadd_pd(zn, c23, temp);
    temp2 = _mm256_mul_pd(temp, temp2);
    temp2 = _mm256_add_pd(temp2, temp2);
    __m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
    __m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
    w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

    return w;
}

// ========== Near Branch ==========
__m256d AddEm(__m256d x)
{
    __m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
    __m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);

    return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

template <size_t MaxOrder>
static __m256d NearBranchSeries(__m256d p)
{
    static constexpr double P[] = {
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

    // Evaluate polynomial using Horner's Method
    __m256d value = _mm256_set1_pd(P[MaxOrder]);
    for (size_t i = 0; i < MaxOrder; i++)
        value = _mm256_fmadd_pd(value, p, _mm256_set1_pd(P[MaxOrder - 1 - i]));

    return value;
}

static __m256d NearBranchW0(__m256d x)
{
    static constexpr double s2e = 2.331643981597124;
    __m256d p = _mm256_mul_pd(_mm256_sqrt_pd(AddEm(x)), _mm256_set1_pd(s2e));

    return NearBranchSeries<24>(p);
}

// ========== Main Function ==========
__m256d MuirpairW0(__m256d x)
{
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.34100), LESS);
    uint32_t nearBranchMask = _mm256_movemask_pd(isNearBranch);

    __m256d result;
    switch (nearBranchMask)
    {
    case 0b0000: [[likely]]
        result = GeneralW0(x);
        break;
    case 0b1111: [[unlikely]]
        result = NearBranchW0(x);
        break;
    default: [[unlikely]]
        result = _mm256_blendv_pd(GeneralW0(x), NearBranchW0(x), isNearBranch);
        break;
    }

    // Fix infinity
    __m256d infinity = _mm256_set1_pd(std::numeric_limits<double>::infinity());
    result = _mm256_blendv_pd(result, infinity, _mm256_cmp_pd(x, infinity, EQUAL));

    // Fix zero
    __m256d isZero = _mm256_cmp_pd(x, _mm256_setzero_pd(), EQUAL);
    result = _mm256_blendv_pd(result, x, isZero);

    return result;
}
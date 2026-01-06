#include <cstdint>
#include <cfloat>
#include <immintrin.h>

#define LESS 0x11

static inline __m256d Epi64ToPd(__m256i x)
{
    __m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    x = _mm256_permutevar8x32_epi32(x, perm);
    return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
}

static inline __m256d LogAccurate(__m256d x)
{
    // === Constants ===
    __m256d ln2 = _mm256_set1_pd(0.69314718055994530942);
    __m256d one = _mm256_set1_pd(1.0);
    __m256d denormThreshold = _mm256_set1_pd(DBL_MIN);
    __m256d denormScale = _mm256_set1_pd(18446744073709551616.0);
    __m256i denormOffset = _mm256_set1_epi64x(64);
    // =================

    // Extract exponent
    __m256d isDenorm = _mm256_cmp_pd(x, denormThreshold, LESS);
    x = _mm256_blendv_pd(x, _mm256_mul_pd(x, denormScale), isDenorm);
    __m256d xd = _mm256_mul_pd(x, _mm256_set1_pd(1.0 / 0.75));
    __m256i dpunn = _mm256_castpd_si256(xd);
    __m256i exp = _mm256_sub_epi64(_mm256_srli_epi64(dpunn, 52), _mm256_set1_epi64x(1023));

    // Extract mantissa
    __m256i punn = _mm256_castpd_si256(x);
    punn = _mm256_sub_epi64(punn, _mm256_slli_epi64(exp, 52));
    __m256d mantissa = _mm256_castsi256_pd(punn);

    // Fix for denormalized exponent
    exp = _mm256_sub_epi64(exp, _mm256_and_si256(_mm256_castpd_si256(isDenorm), denormOffset));

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

static inline __m256d AddEm(__m256d x)
{
	__m256d emHigh = _mm256_set1_pd(0.36787944117144232160);
	__m256d emLow = _mm256_set1_pd(-1.2428753672788363168e-17);

	return _mm256_add_pd(_mm256_add_pd(x, emHigh), emLow);
}

// [EM, -0.277879441171]
static inline __m256d NearBranchWm1(__m256d x)
{
	static constexpr double P[] = {
        -1.000000000000000184466216,-2.331643981596677769194561,-1.812187885818266586169491,-1.936631086156139413432747,-2.353553564102234010397186,-3.06673924788727042848621,-4.179344454070464099829718,-5.764202408685384529372391,-9.991799959165544566758904,7.772856089469904496805629,-208.143553505030026629325,1344.145189003117866322081,-7571.576057874906736835845,31232.91784598825033270358,-96986.89700202980067185836,216959.2797827161006792247,-334089.1543888774939916275,316915.0261048373092021537,-142989.3312531372659257569
	};

	__m256d t = _mm256_sqrt_pd(AddEm(x));

	__m256d numer = _mm256_set1_pd(P[18]);
	for (size_t i = 0; i < 18; i++)
		numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[17 - i]));

    return numer;
}

// [-0.277879441171, -4.1399377188e-8]
static inline __m256d Approx1(__m256d t, __m256d lx)
{
	static constexpr double P[] = {
        8235.87902260375964,-6796.86008640717252,-12443.46895966978809,-6755.27589574189008,-1712.515332068415696,-203.4276623868737738,-8.30812032921245585,-0.02889732975016652775,0.0001706041325417818548
	};
	static constexpr double Q[] = {
        12365.894604528376543058,16052.409820432989815787,9213.150127857532848386,2845.0365311487990888792,477.35031723960426566398,38.336831703708909881334,1.
	};

	__m256d numer = _mm256_set1_pd(P[8]);
	for (size_t i = 0; i < 8; i++)
		numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[7 - i]));

    __m256d denom = _mm256_set1_pd(Q[6]);
    for (size_t i = 0; i < 6; i++)
        denom = _mm256_fmadd_pd(denom, t, _mm256_set1_pd(Q[5 - i]));

    __m256d Y = _mm256_set1_pd(0.666015625);
	return _mm256_add_pd(_mm256_sub_pd(_mm256_div_pd(numer, denom), Y), lx);
}

// [-4.1399377188e-8, 0]
static inline __m256d Approx2(__m256d t, __m256d lx)
{
    static constexpr double P[] = {
        1.5345017051966636633904e6,818003.3778741981291433,79545.29566105822257325,-46046.38588902228941374,-9920.36315181685141171,-549.716885766328377834,-8.1390767501577030303,-0.010035795798761038978,0.00002549847367232381794,-5.0490220466432423419067e-8
    };
    static constexpr double Q[] = {
        512345.13357071709947121,514996.99913940413156071,212748.8328480861983791,41167.296281339095907517,3446.7539760424123748652,110.08635662090830453464,1.
    };

    __m256d numer = _mm256_set1_pd(P[9]);
    for (size_t i = 0; i < 9; i++)
        numer = _mm256_fmadd_pd(numer, t, _mm256_set1_pd(P[8 - i]));

    __m256d denom = _mm256_set1_pd(Q[6]);
    for (size_t i = 0; i < 6; i++)
        denom = _mm256_fmadd_pd(denom, t, _mm256_set1_pd(Q[5 - i]));

    __m256d Y = _mm256_set1_pd(2.9951171875);
    return _mm256_add_pd(_mm256_sub_pd(_mm256_div_pd(numer, denom), Y), lx);
}

static inline __m256d GeneralWm1(__m256d x)
{
    __m256d useFirst = _mm256_cmp_pd(x, _mm256_set1_pd(-4.1399377188e-8), LESS);
    uint32_t useFirstMask = _mm256_movemask_pd(useFirst);

    __m256d negOne = _mm256_set1_pd(-1.0);
    __m256d negX = _mm256_sub_pd(_mm256_setzero_pd(), x);
    __m256d lx = LogAccurate(negX);
    __m256d t = _mm256_sqrt_pd(_mm256_sub_pd(negOne, lx));

    __m256d result;
    switch (useFirstMask)
    {
    case 0b0000:
        result = Approx2(t, lx);
        break;
    case 0b1111:
        result = Approx1(t, lx);
        break;
    default:
        result = _mm256_blendv_pd(Approx2(t, lx), Approx1(t, lx), useFirst);
        break;
    }

    return result;
}

__m256d MuirWm1(__m256d x)
{
    __m256d isNearBranch = _mm256_cmp_pd(x, _mm256_set1_pd(-0.277879441171), LESS);
    uint32_t nearBranchMask = _mm256_movemask_pd(isNearBranch);

    __m256d result;
    switch (nearBranchMask)
    {
    case 0b0000: [[likely]]
        result = GeneralWm1(x);
        break;
    case 0b1111: [[unlikely]]
        result = NearBranchWm1(x);
        break;
    default: [[unlikely]]
        result = _mm256_blendv_pd(GeneralWm1(x), NearBranchWm1(x), isNearBranch);
        break;
    }

    return result;
}
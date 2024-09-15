#include <random>
#include <format>

#include "ulp.h"
#include "wrappers.h"
#include "ReciprocalDistributionEx.h"
#include "histogram.h"

#include "MuirLambertW.h"
#include "reference/ReferenceLambertW.h"
#include "../bench/others/BarryLambertW.h"
#include "boost/math/special_functions/lambert_w.hpp"

#define GREATER 0x1E

double BoostLambertW0(double x)
{
	return boost::math::lambert_w0(x);
}

static inline __m256d Epi64ToPd(__m256i x)
{
	__m256i perm = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
	x = _mm256_permutevar8x32_epi32(x, perm);
	return _mm256_cvtepi32_pd(_mm256_castsi256_si128(x));
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

__m256d Test(__m256d x)
{
	auto [inf, sup] = ReferenceLambertWM1(x[0]);

	__m256d w = _mm256_set1_pd(inf * (1 + pow(10, -7.77534673745)));

	// === Constants ===
	__m256d c23 = _mm256_set1_pd(2.0 / 3.0);
	__m256d one = _mm256_set1_pd(1.0);
	__m256d smallThreshold = _mm256_set1_pd(-1e-300);
	__m256d smallScale = _mm256_set1_pd(4611686018427387904.0);     // 2^62
	__m256d smallOffset = _mm256_set1_pd(42.975125194716609184);    // ln(2^62)
	// =================

	// === Fritsch Iteration ===
	// Compute zn
	__m256d isSmall = _mm256_cmp_pd(x, smallThreshold, GREATER);
	__m256d xScale = _mm256_blendv_pd(x, _mm256_mul_pd(x, smallScale), isSmall);
	__m256d xow = _mm256_div_pd(xScale, w);
	__m256d zn = LogAccurate(xow);
	zn = _mm256_blendv_pd(zn, _mm256_sub_pd(zn, smallOffset), isSmall);
	zn = _mm256_sub_pd(zn, w);

	__m256d temp = _mm256_add_pd(w, one);
	__m256d temp2 = _mm256_fmadd_pd(zn, c23, temp);
	temp2 = _mm256_mul_pd(temp, temp2);
	temp2 = _mm256_add_pd(temp2, temp2);
	__m256d temp3 = _mm256_sub_pd(temp2, _mm256_add_pd(zn, zn));
	__m256d temp4 = _mm256_mul_pd(_mm256_div_pd(zn, temp), _mm256_sub_pd(temp2, zn));
	w = _mm256_mul_pd(w, _mm256_add_pd(_mm256_div_pd(temp4, temp3), one));

	return w;
}

int main()
{
	/*
	static std::mt19937_64 gen{ std::random_device{}() };
	std::uniform_real_distribution<double> dist{ EM_UP, -0.25 };

	double worst = -1;
	for (;;)
	{
		double x = dist(gen);
		auto [inf, sup] = ReferenceLambertWM1(x);

		double approx = MakeSerial<MuirpairWm1>(x);

		uint64_t err = std::max(ULPDistance(approx, inf), ULPDistance(approx, sup));

		if (err <= 4)
			continue;

		if (x > worst)
		{
			worst = x;
			dist = std::uniform_real_distribution<double>{ worst, -0.25 };
			std::cout << std::format("{}\n", worst);
		}
	}
	*/

	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx dist{ EM_UP, 0, false};

	MaxULPRounded(ReferenceLambertWM1, MakeSerial<MuirpairWm1>, [&]() { return dist(gen); });
}
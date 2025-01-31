#include "MuirFukushima.h"

#include <cmath>
#include <limits>
#include <numbers>

#include "MuirFukushimaConstants.h"

template <typename FloatTy>
static inline FloatTy SchroderStep(FloatTy w, FloatTy y)
{
	FloatTy f0 = w - y;
	FloatTy f1 = 1 + y;
	FloatTy f00 = f0 * f0;
	FloatTy f11 = f1 * f1;
	FloatTy f0y = f0 * y;
	return w - 4 * f0 * (6 * f1 * (f11 + f0y) + f00 * y) /
		(f11 * (24 * f11 + 36 * f0y) +
			f00 * (6 * y * y + 8 * f1 * y + f0y));
}

template <typename FloatTy, size_t Order>
static inline FloatTy NearBranchSeries(FloatTy p)
{
	static constexpr FloatTy P[] = {
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
		-0.0000581136075044138168
	};

	FloatTy numer = P[Order];
	for (size_t i = 0; i < Order; i++)
		numer = numer * p + P[Order - 1 - i];

	return numer;
}

template <typename FloatTy>
static inline FloatTy NearZeroSeries(FloatTy x)
{
	static constexpr FloatTy P[] = {
		0,
		1,
		-1,
		1.5,
		-2.6666666666666666667,
		5.2083333333333333333,
		-10.8,
		23.343055555555555556,
		-52.012698412698412698,
		118.62522321428571429,
		-275.57319223985890653,
		649.78717234347442681,
		-1551.1605194805194805,
		3741.4497029592385495,
		-9104.5002411580189358,
		22324.308512706601434,
		-55103.621972903835338,
		136808.86090394293563
	};

	constexpr size_t Order = (std::is_same_v<FloatTy, double>) ? 17 : 7;

	FloatTy numer = P[Order];
	for (size_t i = 0; i < Order; i++)
		numer = numer * x + P[Order - 1 - i];

	return numer;
}

template <typename FloatTy>
static inline FloatTy NearBranchW0(FloatTy x)
{
	constexpr FloatTy e2 = std::numbers::e_v<FloatTy> * 2;
	FloatTy p = std::sqrt(e2 * x + 2);
	if constexpr (std::is_same_v<FloatTy, double>)
	{
		if (std::abs(p) < 0.01159) return NearBranchSeries<FloatTy, 6>(p);
		if (std::abs(p) > 0.0766) return NearBranchSeries<FloatTy, 20>(p);
	}
	return NearBranchSeries<FloatTy, 10>(p);
}

float MuirFukushimaW0(float x)
{
	return 0.0f;
}

int IntegerPartW0(double x)
{
	DECLARE_W0_G;

	// Initial linear search for small x
	int n;
	for (n = 0; n <= 2; ++n)
		if (G[n] > x)
			return n - 1;

	// Search at powers of 2
	n = 2;
	for (int j = 1; j <= 5; ++j)
	{
		n *= 2;
		if (G[n] > x) break;
	}

	// Bisect to determine final n
	int nh = n / 2;
	for (int j = 1; j <= 5; ++j)
	{
		nh /= 2;
		if (nh <= 0)
			break;
		if (G[n - nh] > x)
			n -= nh;
	}

	return n - 1;
}

double MuirFukushimaW0(double x)
{
	DECLARE_W0_E;
	DECLARE_W0_A;
	DECLARE_W0_B;

	// Edge cases
	if (abs(x) < 0.05) return NearZeroSeries(x);
	if (x < -0.35) return NearBranchW0(x);

	// Get integer part
	int n = IntegerPartW0(x);

	// Determine number of bisections
	int jmax = 8;
	if (x <= -0.36) jmax = 12;
	else if (x <= -0.3) jmax = 11;
	else if (n <= 0) jmax = 10;
	else if (n <= 1) jmax = 9;

	// Bisection
	double y = x * E[n + 1];
	double w = n;
	for (int j = 0; j < jmax; ++j)
	{
		const double wj = w + B[j];
		const double yj = y * A[j];
		if (wj < yj)
		{
			w = wj;
			y = yj;
		}
	}

	return SchroderStep(w, y);
}

float MuirFukushimaWm1(float x)
{
	return 0.0f;
}

double MuirFukushimaWm1(double x)
{
	return 0.0;
}
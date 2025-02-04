#include "MuirFukushima.h"

#include <cmath>
#include <limits>
#include <numbers>

#include <immintrin.h>

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

template <typename FloatTy, bool IsWm1>
static inline FloatTy NearBranch(FloatTy x)
{
	constexpr FloatTy e2 = std::numbers::e_v<FloatTy> * 2;
	FloatTy p = std::sqrt(e2 * x + 2);
	if constexpr (IsWm1) p = -p;
	if constexpr (std::is_same_v<FloatTy, double>)
	{
		if (std::abs(p) < 0.01159) return NearBranchSeries<FloatTy, 6>(p);
		if (std::abs(p) > 0.0766) return NearBranchSeries<FloatTy, 20>(p);
	}
	return NearBranchSeries<FloatTy, 10>(p);
}

// ===== W0f =====
static inline int IntegerPartW0f(float x)
{
	DECLARE_W0F_G;

	// Initial linear search for small x
	int n;
	for (n = 0; n <= 2; ++n)
		if (G[n] > x)
			return n - 1;

	// Search at powers of 2
	n = 2;
	for (int j = 1; j <= 4; ++j)
	{
		n *= 2;
		if (G[n] > x) break;
	}

	// Bisect to determine final n
	int nh = n / 2;
	for (int j = 1; j <= 4; ++j)
	{
		nh /= 2;
		if (nh <= 0) break;
		if (G[n - nh] > x) n -= nh;
	}

	return n - 1;
}

template <size_t StartIter, size_t NumIter>
static inline void BisectionW0f(float* wp, float* yp, float wi, float yi)
{
	DECLARE_W0F_A;

	__m128 y = _mm_set_ss(yi);
	__m128 w = _mm_set_ss(wi);
	__m128 half = _mm_set_ss(0.5);
	constexpr float bStart = 0.5 / (1 << StartIter);
	__m128 b = _mm_set_ss(bStart);
	for (size_t j = 0; j < NumIter; j++)
	{
		__m128 wj = _mm_add_ss(w, b);
		__m128 yj = _mm_mul_ss(y, _mm_set_ss(A[StartIter + j]));
		w = _mm_blendv_ps(w, wj, _mm_cmplt_ss(wj, yj));
		y = _mm_blendv_ps(y, yj, _mm_cmplt_ss(wj, yj));
		b = _mm_mul_ss(b, half);
	}

	*wp = w[0];
	*yp = y[0];
}

float MuirFukushimaW0(float x)
{
	DECLARE_W0F_E;

	// Edge cases
	if (x >= 2.5268146e+15f) return std::numeric_limits<float>::quiet_NaN();
	if (std::abs(x) < 0.05f) return NearZeroSeries(x);
	if (x < -0.33f) return NearBranch<float, false>(x);

	// Get integer part
	int n = IntegerPartW0f(x);

	// Do the initial 2 bisections
	float w = n;
	float y = x * E[n + 1];
	BisectionW0f<0, 2>(&w, &y, w, y);

	// Perform remaining bisections if necessary
	if (x <= -0.3f) BisectionW0f<2, 3>(&w, &y, w, y);
	else if (x <= -0.1f) BisectionW0f<2, 2>(&w, &y, w, y);
	else if (n < 8) BisectionW0f<2, 1>(&w, &y, w, y);

	return SchroderStep(w, y);
}

// ===== W0 =====
static inline int IntegerPartW0(double x)
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
		if (nh <= 0) break;
		if (G[n - nh] > x) n -= nh;
	}

	return n - 1;
}

template <size_t StartIter, size_t NumIter>
static inline void BisectionW0(double* wp, double* yp, double wi, double yi)
{
	DECLARE_W0_A;

	__m128d y = _mm_set_sd(yi);
	__m128d w = _mm_set_sd(wi);
	__m128d half = _mm_set_sd(0.5);
	constexpr double bStart = 0.5 / (1 << StartIter);
	__m128d b = _mm_set_sd(bStart);
	for (size_t j = 0; j < NumIter; j++)
	{
		__m128d wj = _mm_add_sd(w, b);
		__m128d yj = _mm_mul_sd(y, _mm_set_sd(A[StartIter + j]));
		w = _mm_blendv_pd(w, wj, _mm_cmplt_sd(wj, yj));
		y = _mm_blendv_pd(y, yj, _mm_cmplt_sd(wj, yj));
		b = _mm_mul_sd(b, half);
	}

	*wp = w[0];
	*yp = y[0];
}

double MuirFukushimaW0(double x)
{
	DECLARE_W0_E;

	// Edge cases
	if (x >= 3.990495411719435e+29) return std::numeric_limits<double>::quiet_NaN();
	if (abs(x) < 0.05) return NearZeroSeries(x);
	if (x < -0.35) return NearBranch<double, false>(x);

	// Get integer part
	int n = IntegerPartW0(x);

	// Do the initial 8 bisections
	double w = n;
	double y = x * E[n + 1];
	BisectionW0<0, 8>(&w, &y, w, y);

	// Perform remaining bisections if necessary
	if (x <= -0.3) BisectionW0<8, 3>(&w, &y, w, y);
	else if (n <= 0) BisectionW0<8, 2>(&w, &y, w, y);
	else if (n <= 1) BisectionW0<8, 1>(&w, &y, w, y);

	return SchroderStep(w, y);
}

// ===== Wm1f =====
static inline int IntegerPartWm1f(float x)
{
	DECLARE_WM1F_G;

	if (G[1] > x) return 1;

	// Search at powers of 2
	int n = 2;
	for (int j = 1; j <= 4; ++j)
	{
		n *= 2;
		if (G[n - 1] > x) break;
	}

	// Bisect to determine final n
	int nh = n / 2;
	for (int j = 1; j <= 4; ++j)
	{
		nh /= 2;
		if (nh <= 0) break;
		if (G[n - nh - 1] > x) n -= nh;
	}

	return n - 1;
}

template <size_t StartIter, size_t NumIter>
static inline void BisectionWm1f(float* wp, float* yp, float wi, float yi)
{
	DECLARE_WM1F_A;

	__m128 y = _mm_set_ss(yi);
	__m128 w = _mm_set_ss(wi);
	__m128 half = _mm_set_ss(0.5);
	constexpr float bStart = 0.5 / (1 << StartIter);
	__m128 b = _mm_set_ss(bStart);
	for (size_t j = 0; j < NumIter; j++)
	{
		__m128 wj = _mm_sub_ss(w, b);
		__m128 yj = _mm_mul_ss(y, _mm_set_ss(A[StartIter + j]));
		w = _mm_blendv_ps(w, wj, _mm_cmplt_ss(wj, yj));
		y = _mm_blendv_ps(y, yj, _mm_cmplt_ss(wj, yj));
		b = _mm_mul_ss(b, half);
	}

	*wp = w[0];
	*yp = y[0];
}

float MuirFukushimaWm1(float x)
{
	DECLARE_WM1F_E;

	// Edge cases
	if (x >= -4.052533e-13f) return std::numeric_limits<float>::quiet_NaN();
	if (x < -0.33f) return NearBranch<float, true>(x);

	// Get integer part
	int n = IntegerPartWm1f(x);

	// Do the initial 2 bisections
	float w = -n;
	float y = x * E[n - 1];
	BisectionWm1f<0, 2>(&w, &y, w, y);

	// Perform remaining bisections if necessary
	if (n <= 1) BisectionWm1f<2, 3>(&w, &y, w, y);
	else if (n <= 2) BisectionWm1f<2, 2>(&w, &y, w, y);
	else if (n <= 7) BisectionWm1f<2, 1>(&w, &y, w, y);

	return SchroderStep(w, y);
}

// ===== Wm1 =====
static inline int IntegerPartWm1(double x)
{
	DECLARE_WM1_G;

	if (G[1] > x) return 1;

	// Search at powers of 2
	int n = 2;
	for (int j = 1; j <= 5; ++j)
	{
		n *= 2;
		if (G[n - 1] > x) break;
	}

	// Bisect to determine final n
	int nh = n / 2;
	for (int j = 1; j <= 5; ++j)
	{
		nh /= 2;
		if (nh <= 0) break;
		if (G[n - nh - 1] > x) n -= nh;
	}

	return n - 1;
}

template <size_t StartIter, size_t NumIter>
static inline void BisectionWm1(double* wp, double* yp, double wi, double yi)
{
	DECLARE_WM1_A;

	__m128d y = _mm_set_sd(yi);
	__m128d w = _mm_set_sd(wi);
	__m128d half = _mm_set_sd(0.5);
	constexpr double bStart = 0.5 / (1 << StartIter);
	__m128d b = _mm_set_sd(bStart);
	for (size_t j = 0; j < NumIter; j++)
	{
		__m128d wj = _mm_sub_sd(w, b);
		__m128d yj = _mm_mul_sd(y, _mm_set_sd(A[StartIter + j]));
		w = _mm_blendv_pd(w, wj, _mm_cmplt_sd(wj, yj));
		y = _mm_blendv_pd(y, yj, _mm_cmplt_sd(wj, yj));
		b = _mm_mul_sd(b, half);
	}

	*wp = w[0];
	*yp = y[0];
}

double MuirFukushimaWm1(double x)
{
	DECLARE_WM1_E;

	// Edge cases
	if (x >= -1.0264389699511283e-26) return std::numeric_limits<double>::quiet_NaN();
	if (x < -0.35) return NearBranch<double, true>(x);

	// Get integer part
	int n = IntegerPartWm1(x);

	// Do the initial 8 bisections
	double w = -n;
	double y = x * E[n - 1];
	BisectionWm1<0, 8>(&w, &y, w, y);

	// Perform remaining bisections if necessary
	if (n <= 1) BisectionWm1<8, 3>(&w, &y, w, y);
	else if (n <= 2) BisectionWm1<8, 2>(&w, &y, w, y);
	else if (n <= 7) BisectionWm1<8, 1>(&w, &y, w, y);

	return SchroderStep(w, y);
}
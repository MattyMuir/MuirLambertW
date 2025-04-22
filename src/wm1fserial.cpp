#include <cmath>

static inline float NearBranchWm1(float x)
{
	static constexpr double e2 = 5.43656365691809;
	double p = sqrt(e2 * x + 2.0);

	// Polynomial approximation coefficients for algorithm 7, index 1, order 5
	static constexpr double P[] = {
		-0.9999999811456308,
		-1.0000044654258957,
		-0.3331646630001799,
		-0.1550676095217873,
		-0.06595497662216744,
		-0.07987153028272707
	};

	double res = P[5];
	for (size_t i = 0; i < 5; i++)
		res = res * p + P[4 - i];

	return res;
}

static inline float GeneralWm1(float x)
{
	// Rational approximation coefficients for algorithm 7, index 2, order 5/3
	static constexpr double P[] = {
		-87.77561192149614,
		-184.9753116484922,
		-159.16694306575346,
		-69.07056916463522,
		-14.653530143154867,
		-0.999363416489587
	};
	static constexpr double Q[] = {
		87.77537609045154,
		60.84479581959989,
		14.59140282955174,
		1
	};

	double t = sqrt(-1 - log((double)-x));

	double numer = P[5];
	for (size_t i = 0; i < 5; i++)
		numer = numer * t + P[4 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * t + Q[2 - i];

	return numer / denom;
}

float MuirWm1(float x)
{
	return (x < -0.3498715f) ? NearBranchWm1(x) : GeneralWm1(x);
}
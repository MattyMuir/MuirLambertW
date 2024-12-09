#include <cmath>

static inline float NearBranchWm1(float x)
{
	static constexpr double e2 = 5.43656365691809;
	double p = sqrt(e2 * x + 2.0);

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
	static constexpr double P[] = {
		-246.9023719292306,
		-368.15952949040997,
		-224.1791808874142,
		-68.85182301857488,
		-10.34301280940623,
		-0.4996779961653583
	};

	static constexpr double Q[] = {
		246.90158428606017,
		121.26384144163241,
		20.59776676808664,
		1
	};

	double t = sqrt(-2 - 2 * log((double)-x));

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
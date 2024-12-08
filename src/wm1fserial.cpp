#include <cmath>

static inline float NearBranchWm1(float x)
{
	static constexpr double s2e = 2.331643981597124;

	double p = s2e * sqrt((double)x + 0.36787944117144233);

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
		-2101.555169658076,
		-3413.0457024602106,
		-2345.4071921263444,
		-864.1804177336671,
		-175.99964384176346,
		-17.64071303855079,
		-0.4998769261313046
	};

	static constexpr double Q[] = {
		2101.5551872949245,
		1311.4898275251383,
		333.4030604186147,
		35.228646667156625,
		1
	};

	double t = sqrt(-2 - 2 * log((double)-x));

	double numer = P[6];
	for (size_t i = 0; i < 6; i++)
		numer = numer * t + P[5 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

float MuirWm1(float x)
{
	return (x < -0.3498715f) ? NearBranchWm1(x) : GeneralWm1(x);
}
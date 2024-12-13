#include <cmath>

static inline float FirstApproxW0(float x)
{
	static constexpr double P[] = {
		0,
		152.91909285006932,
		1024.8955779609623,
		2455.092024342446,
		2528.414753137949,
		1066.824695440063,
		148.54000057530493,
		4.0539143590175835
	};

	static constexpr double Q[] = {
		152.91909349763372,
		1177.8147226529784,
		3403.5283401899756,
		4573.002542975176,
		2878.913736848507,
		761.147965396105,
		65.5632478027466,
		1
	};

	double numer = P[7];
	for (size_t i = 0; i < 7; i++)
		numer = numer * x + P[6 - i];

	double denom = Q[7];
	for (size_t i = 0; i < 7; i++)
		denom = denom * x + Q[6 - i];

	return numer / denom;
}

static inline float SecondApproxW0(float x)
{
	static constexpr double P[] = {
		245182.20097823755,
		280243.5212428723,
		142843.813324628,
		40353.72076097795,
		5776.914448840662,
		184.83613670644033,
		0.9984483567344636
	};

	static constexpr double Q[] = {
		432788.26007218857,
		216948.13159273885,
		58081.26591912717,
		6594.751582203545,
		191.21022696372594,
		1
	};

	double t = log((double)x);

	double numer = P[6];
	for (size_t i = 0; i < 6; i++)
		numer = numer * t + P[5 - i];

	double denom = Q[5];
	for (size_t i = 0; i < 5; i++)
		denom = denom * t + Q[4 - i];

	return numer / denom;
}

static inline float NearBranchW0(float x)
{
	static constexpr double e2 = 5.43656365691809;

	static constexpr double P[] = {
		-0.9999999781289544,
		0.9999966080647236,
		-0.33324531164727067,
		0.15189891604646868,
		-0.07530393941472714,
		0.03290035332102544,
		-0.008369773627101843
	};

	double p = sqrt(e2 * x + 2.0);

	double res = P[6];
	for (size_t i = 0; i < 6; i++)
		res = res * p + P[5 - i];

	return res;
}

float MuirW0(float x)
{
	//return SecondApproxW0(x);
	return (x < -0.3f) ? NearBranchW0(x) : ((x < 6.9035267829895019531f) ? FirstApproxW0(x) : SecondApproxW0(x));
}
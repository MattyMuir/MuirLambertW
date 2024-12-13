#include <cmath>

static inline float FirstApproxW0(float x)
{
	static constexpr double P[] = {
		0,
		165.51561672164559,
		1104.9153130867758,
		2632.284078577963,
		2689.464120405435,
		1121.2923665114324,
		153.3374641092571,
		4.077322829553558
	};

	static constexpr double Q[] = {
		165.51561558818844,
		1270.4310030077481,
		3654.442208397931,
		4879.631928655197,
		3045.0058891120098,
		794.8712729472717,
		67.22857835896016,
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

	double t = std::logf(x);

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
	return (x < -0.3f) ? NearBranchW0(x) : ((x < 7.38905609893f) ? FirstApproxW0(x) : SecondApproxW0(x));
}
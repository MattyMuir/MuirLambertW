#include <cmath>

// [EM, -0.00000520258255623]
static inline float Approx1(float x)
{
	static constexpr double P[] = {
		-0.999999981377920024398681,0.08222516843865509649862124,0.3767400006110998370390504,0.04604539061352127522237554
	};
	static constexpr double Q[] = {
		1,0.9177734753926605706900684,0.2077166606126629054461969,0.008444149579041735852782045
	};

	static constexpr double e2 = 5.43656365691809;
	double p = sqrt(e2 * x + 2.0);

	double numer = P[3];
	for (size_t i = 0; i < 3; i++)
		numer = numer * p + P[2 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * p + Q[2 - i];

	return numer / denom;
}

// [-0.00000520258255623, 0.9]
static inline float Approx2(float x)
{
	static constexpr double P[] = {
		-0.9999999678909199544267152,-2.584348453493142588704776,-1.294389423480461386526167,-0.01005043244141658450465425
	};
	static constexpr double Q[] = {
		1,4.084341700580141865622986,4.754463651441859052466608,1.455430947893811785272415
	};

	double numer = P[3];
	for (size_t i = 0; i < 3; i++)
		numer = numer * x + P[2 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * x + Q[2 - i];

	return x + (numer / denom) * x * x;
}

// [0.9, 8.633]
static inline float Approx3(float x)
{
	static constexpr double P[] = {
		0.001916675355730279319460353,0.9853176179570127833014898,1.098324712410893608271392,0.2144434781956386510863834,0.00606866731378611170881161
	};
	static constexpr double Q[] = {
		1,2.039365418315944583483287,0.9280500438839044664089006,0.09722552755883852378797482,0.00147642461414064698231248
	};

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * x + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * x + Q[3 - i];

	return numer / denom;
}

// [8.633, 72.2888]
static inline float Approx4(float x)
{
	static constexpr double P[] = {
		0.1696617959997037522279642,0.6415324632315349922441918,0.09539854751430477187818722,0.002310981226848268753976592,7.664599465549063222159969e-6
	};
	static constexpr double Q[] = {
		1,0.4975013559377048349066393,0.03989671017187767773553451,0.0006277408159739880464674528,1.322993800985472670212182e-6
	};

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * x + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * x + Q[3 - i];

	return numer / denom;
}

// [72.2888, 5473.26234225]
static inline float Approx5(float x)
{
	static constexpr double P[] = {
		-0.4376440913873589471746691,1.084084375199522288100761,0.168206363802967484939933,0.004014540398865475546010566,0.00001306140558267055547841221
	};
	static constexpr double Q[] = {
		1,0.4428741699016309115976059,0.03324076358342433803352941,0.0005021761267753373851834467,1.033871460909154906803944e-6
	};

	double t = sqrt((double)x);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [5473.26234225, 3545455.51184]
static inline float Approx6(float x)
{
	static constexpr double P[] = {
		1.449957995853132329744377,0.4356423156126395218718242,0.007945369031738161484065351,0.00002811644744126273715215578,2.12912180173491327898596e-8,2.399796274848595480711734e-12
	};
	static constexpr double Q[] = {
		1,0.08029225944656790261737529,0.000963083352324697374868256,2.582092561765425030733291e-6,1.532333581966206535089403e-9,1.287347360469109661707269e-13
	};

	double t = sqrt((double)x);

	double numer = P[5];
	for (size_t i = 0; i < 5; i++)
		numer = numer * t + P[4 - i];

	double denom = Q[5];
	for (size_t i = 0; i < 5; i++)
		denom = denom * t + Q[4 - i];

	return numer / denom;
}

// [3545455.51184, 1.8750526518e10]
static inline float Approx7(float x)
{
	static constexpr double P[] = {
		1.19025368111793336039939,1.342550099584281646786995,0.03991431268170822350826933,0.0001886561640009630930369888,1.197404437845505041826168e-7
	};
	static constexpr double Q[] = {
		1,0.1328548982090831343746002,0.002412434392835541988028041,8.149906457806579395898324e-6,3.621266492575407455308264e-9
	};

	double t = sqrt(sqrt((double)x));

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [1.8750526518e10, FLT_MAX]
static inline float Approx8(float x)
{
	static constexpr double P[] = {
		0.3599923913218089460631521,0.482340264345297233991295,0.1634673946984464208758099,0.002923429305941059114264689
	};
	static constexpr double Q[] = {
		1,0.1791274390047670094491925,0.002937246705565411862651593,-1.798181964681319281352972e-8
	};

	double t = log((double)x);

	double numer = P[3];
	for (size_t i = 0; i < 3; i++)
		numer = numer * t + P[2 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * t + Q[2 - i];

	return numer / denom;
}

float MuirW0v2(float x)
{
	if (x < 72.2888f)
	{
		if (x < 0.9f)
		{
			if (x < -0.00000520258255623f) return Approx1(x);
			return Approx2(x);
		}
		if (x < 8.633f) return Approx3(x);
		return Approx4(x);
	}
	if (x < 3545455.51184f)
	{
		if (x < 5473.26234225f) return Approx5(x);
		return Approx6(x);
	}
	if (x < 1.8750526518e10f) return Approx7(x);
	return Approx8(x);
}
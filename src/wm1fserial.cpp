#include <cmath>

// [EM, -0.103006243528]
static inline float Approx1(float x)
{
	static constexpr double P[] = {
		-1.0000000289164983219173,
		0.7484377282132577749822601,
		0.3718673404409796721564797,
		-0.3760756806217422409838732,
		0.06188486908948103343472952
	};
	static constexpr double Q[] = {
		1,
		-1.748435627144009813001159,
		1.043209266905849599870069,
		-0.2369809021813308821994934,
		0.0145973490111531148723476
	};

	static constexpr double e2 = 5.43656365691809;
	double p = sqrt(e2 * x + 2.0);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * p + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * p + Q[3 - i];

	return numer / denom;
}

// [-0.103006243528, -0.00423944477838]
static inline float Approx2(float x)
{
	static constexpr double P[] = {
		-1.246610243659013302100792,
		-1.840405389569920475830382,
		-5.057542400159065240856384,
		-1.609661227270221046571363,
		-0.05174401895621838370850601
	};
	static constexpr double Q[] = {
		1,
		1.34572414751901002995512,
		1.127467170560772841677794,
		0.1839003550545753464700899,
		0.003563775943212132799644493
	};

	double t = sqrt(-1.0 - 1.0 / (x * 2.718281828459045));

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.00423944477838, -0.0000523750621517]
static inline float Approx3(float x)
{
	static constexpr double P[] = {
		-0.2739876278365825998500854,
		-3.586369508785214811994107,
		-0.4601152729439305196776908,
		-0.009461127078624370169441965,
		-0.00002625200058912045255396949
	};
	static constexpr double Q[] = {
		1,
		0.5742709990602670957507459,
		0.04568777602546817846367936,
		0.0006847238651296436794073281,
		1.362664148332156080899477e-6
	};

	double t = sqrt(-1.0 - 1.0 / (x * 2.718281828459045));

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.0000523750621517, -1.1274536796e-10]
static inline float Approx4(float x)
{
	static constexpr double P[] = {
		0.6460922708494494169644241,
		-5.63606810723108947436321,
		-0.8523673575412450020053934,
		-0.02393514755644661139781376,
		-0.0001433721837039804358904018,
		-1.282465374370997827518384e-7
	};
	static constexpr double Q[] = {
		1,
		0.5746573260180125538160536,
		0.05159856094451173125023817,
		0.001061411130045904009412654,
		4.901386251165684872301854e-6,
		3.238895413874947261010569e-9
	};

	double t = sqrt(1.0 / sqrt(x * -2.718281828459045) - 1.0);

	double numer = P[5];
	for (size_t i = 0; i < 5; i++)
		numer = numer * t + P[4 - i];

	double denom = Q[5];
	for (size_t i = 0; i < 5; i++)
		denom = denom * t + Q[4 - i];

	return numer / denom;
}

// [-1.1274536796e-10, 0]
static inline float Approx5(float x)
{
	static constexpr double P[] = {
		-1.00075040133045730502965,
		-1.696302052219855127272343,
		-1.042001775267366896884873,
		-0.2739406543277726856340339
	};
	static constexpr double Q[] = {
		1,
		0.2750686586148924203848048,
		-0.00001834201536276028988413627
	};

	double t = sqrt(-1 - log((double)-x));

	double numer = P[3];
	for (size_t i = 0; i < 3; i++)
		numer = numer * t + P[2 - i];

	double denom = Q[2];
	for (size_t i = 0; i < 2; i++)
		denom = denom * t + Q[1 - i];

	return numer / denom;
}

float MuirWm1(float x)
{
	if (x < -0.00423944477838f)
	{
		if (x < -0.103006243528f) return Approx1(x);
		return Approx2(x);
	}
	if (x < -0.0000523750621517f) return Approx3(x);
	if (x < -1.1274536796e-10f) return Approx4(x);
	return Approx5(x);
}
#include <cmath>

// [EM, -0.103006243528]
static inline float Approx1(float x)
{
	static constexpr double P[] = {
		-2.3178059085766222762974,4.044780548366448857825,4.6858627990364218546778,-11.0494082385605575260139,4.2394594417272972291861
	};
	static constexpr double Q[] = {
		2.31780584155379354981898,-9.4490672345315887428158,13.1453746102011527767953,-6.9626909366033831452371,1.
	};

	static constexpr double em = 0.36787944117144233;
	double t = sqrt(x + em);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.103006243528, -0.015006002401]
static inline float Approx2(float x)
{
	static constexpr double P[] = {
		531540.53921623203536812,202802.2320397417372886,-79332.737287649801637587,2538.9370756728162401734,-9.6701879924297624377769
	};
	static constexpr double Q[] = {
		-6067.30569382733974164204,-122534.013144607968987478,17651.162749470299619527,-378.516081061524253332244,1.
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.015006002401, -0.001643817602]
static inline float Approx3(float x)
{
	static constexpr double P[] = {
		-1.694225608545212417645263,0.523359981688664248531702,-0.00903188397390566799151936,0.00002547986877248439331854359,-9.586128141156491362787871e-9
	};
	static constexpr double Q[] = {
		1.,-0.1011130175686267502294596,0.001243981563287747095602623,-2.755115242048638238578842e-6,7.92886993609830424549085e-10
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.001643817602, -0.000162503859467]
static inline float Approx4(float x)
{
	static constexpr double P[] = {
		-4.822424349679403110422753,0.07695325040727959504917168,-0.0001197058325103768354721113,3.223210222966952612260154e-8,-1.163789406379266609843764e-12
	};
	static constexpr double Q[] = {
		1.,-0.009778624611220878003528704,0.0000121627692142777222771119,-2.738421411479261857250144e-9,7.991687887805711435842493e-14
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.000162503859467, -0.0000149302756129]
static inline float Approx5(float x)
{
	static constexpr double P[] = {
		-7.607081997788888412906741,0.009668777144551208420170272,-1.346392563619916635162224e-6,3.29866644924982480691253e-11,-1.079283699913608198435982e-16
	};
	static constexpr double Q[] = {
		1.,-0.0009201364635333133892089738,1.083200572324020855326532e-7,-2.302570870384622421169944e-12,6.31658460460245123918948e-18
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.0000149302756129, -0.0000012926744141]
static inline float Approx6(float x)
{
	static constexpr double P[] = {
		-10.31853173596190550719555,0.001076789092131532074985554,-1.290405073045279969130336e-8,2.736193034125162726605937e-14,-7.703023769655454811162045e-21
	};
	static constexpr double Q[] = {
		1.,-0.00008190808297505144935720515,8.580378296829680456558841e-10,-1.617601227399575152036968e-15,3.919389404631643161615941e-22
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-0.0000012926744141, -1.0640448601e-7]
static inline float Approx7(float x)
{
	static constexpr double P[] = {
		-13.0150750293473491889842,0.0001093812722434170206011437,-1.081341146001545346565552e-10,1.894727925015306199173538e-17,-4.381899307863175949101766e-25
	};
	static constexpr double Q[] = {
		1.,-6.923162516044875660327393e-6,6.118454921919997483957043e-12,-9.698686550263043307039831e-19,1.968647414006654711449132e-26
	};

	double t = 1.0 / (double)x;

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	return numer / denom;
}

// [-1.0640448601e-7, 0]
static inline float Approx8(float x)
{
	static constexpr double P[] = {
		5114.3511421762101261432,-26247.164076947694742457,5539.2828488128405071193,-190.702382635947266295553,1.00144049343340160162538
	};
	static constexpr double Q[] = {
		-14547.4094153324000071882,4761.76106690377968183025,-184.293152826205802261266,1.
	};

	double t = log(-(double)x);

	double numer = P[4];
	for (size_t i = 0; i < 4; i++)
		numer = numer * t + P[3 - i];

	double denom = Q[3];
	for (size_t i = 0; i < 3; i++)
		denom = denom * t + Q[2 - i];

	return numer / denom;
}

float MuirWm1(float x)
{
	if (x < -0.000162503859467f)
	{
		if (x < -0.015006002401f)
		{
			if (x < -0.103006243528f) return Approx1(x);
			return Approx2(x);
		}
		if (x < -0.001643817602f) return Approx3(x);
		return Approx4(x);
	}
	if (x < -0.0000012926744141f)
	{
		if (x < -0.0000149302756129f) return Approx5(x);
		return Approx6(x);
	}
	if (x < -1.0640448601e-7f) return Approx7(x);
	return Approx8(x);
}
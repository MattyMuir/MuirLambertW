#include <cmath>

static inline double Approx1(double t)
{
	static constexpr double P[] = {
		-0.9999999941229416,
		-1.414213666694153,
		-0.6666658409084293,
		-0.07857125226020724,
		0.01482632573225921,
		-0.0013326922503326796,
		-0.000438877944874127,
		0.0002620093772182476,
		-6.700663061738577e-05,
		7.775177337397004e-06
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx2(double t)
{
	static constexpr double P[] = {
		-0.9999991052138872,
		-1.4142228432466455,
		-0.6666233526751946,
		-0.0786872425299455,
		0.015032421058479434,
		-0.0015802818965254203,
		-0.00023746305896627065,
		0.00015487808233076116,
		-3.31687051339452e-05,
		2.937564905389504e-06
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx3(double t)
{
	static constexpr double P[] = {
		-0.9999837643007623,
		-1.414330676006523,
		-0.6662838511310054,
		-0.07931600956592565,
		0.015787732414707312,
		-0.002190848517035432,
		9.47697556114369e-05,
		3.751635591470898e-05,
		-8.74657331151851e-06,
		6.570100158930688e-07
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx4(double t)
{
	static constexpr double P[] = {
		-0.9999000278350678,
		-1.4147589111597967,
		-0.6653039899082682,
		-0.08063290032438235,
		0.016933488606737712,
		-0.002860128522528426,
		0.0003572412373170384,
		-2.9113773771046054e-05,
		1.1866993772296217e-06,
		-5.412129857849659e-09
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx5(double t)
{
	static constexpr double P[] = {
		-0.999850714771072,
		-1.4149825109462812,
		-0.6648633657862012,
		-0.08113151291240682,
		0.017292127123492998,
		-0.003030678014864718,
		0.0004109783316108222,
		-3.994840432530897e-05,
		2.456628039334284e-06,
		-7.139831850730499e-08
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx6(double t)
{
	static constexpr double P[] = {
		-1.0012901826225655,
		-1.4110383515930187,
		-0.6696841902249538,
		-0.07768098840432018,
		0.0156980103566251,
		-0.0025376305645794985,
		0.0003088729349054402,
		-2.6294660971587267e-05,
		1.3867535014945948e-06,
		-3.396906279967354e-08
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx7(double t)
{
	static constexpr double P[] = {
		-1.0078135058938402,
		-1.3967768836447279,
		-0.6835888710223645,
		-0.06974512555393052,
		0.01277595429983442,
		-0.0018177561703198663,
		0.00019021206252622816,
		-1.3675052534933088e-05,
		6.010411325163889e-07,
		-1.2149730369168439e-08
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx8(double t)
{
	static constexpr double P[] = {
		-1.025493341416611,
		-1.3656030356706288,
		-0.7080945528745963,
		-0.058472491861552074,
		0.009431883158797477,
		-0.0011542948304495465,
		0.00010217967479738935,
		-6.142344137744855e-06,
		2.2388071182088044e-07,
		-3.731062577318834e-09
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx9(double t)
{
	static constexpr double P[] = {
		-1.0623301828314509,
		-1.3127800987007474,
		-0.7418580356029199,
		-0.04584647567539287,
		0.00638762316512179,
		-0.0006635164088702232,
		4.9278069228776416e-05,
		-2.4659439912335426e-06,
		7.44194226154773e-08,
		-1.0229854607093058e-09
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx10(double t)
{
	static constexpr double P[] = {
		-1.1268311651757204,
		-1.2371930463088006,
		-0.7813372548048173,
		-0.03378404606263004,
		0.004011638685185863,
		-0.0003506325737364509,
		2.1733153004770373e-05,
		-9.027717599831511e-07,
		2.253238962115171e-08,
		-2.5549922925469494e-10
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx11(double t)
{
	static constexpr double P[] = {
		-1.2268394843789134,
		-1.141209231982697,
		-0.822394262519135,
		-0.023510841259809026,
		0.0023545395169515884,
		-0.0001719425265225157,
		8.852225596550371e-06,
		-3.0425701532633937e-07,
		6.267120938570492e-09,
		-5.853929332989299e-11
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx12(double t)
{
	static constexpr double P[] = {
		-1.3684016349998525,
		-1.0299502440550017,
		-0.8613670517877131,
		-0.015524906734856384,
		0.0012996187836371727,
		-7.87831578320482e-05,
		3.3526537683325287e-06,
		-9.49853715781347e-08,
		1.6096905444416979e-09,
		-1.2353735397314039e-11
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx13(double t)
{
	static constexpr double P[] = {
		-1.555380238870257,
		-0.909839477001033,
		-0.8957578664454909,
		-0.009764262966383427,
		0.0006775228248998383,
		-3.386913807490711e-05,
		1.1848245564772124e-06,
		-2.7537600442056433e-08,
		3.8230200805706366e-10,
		-2.4011711339170306e-12
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx14(double t)
{
	static constexpr double P[] = {
		-1.7899059270488393,
		-0.7871624743141583,
		-0.9243642132996572,
		-0.005861396119276726,
		0.0003341906889100253,
		-1.3674689384291126e-05,
		3.906546192497785e-07,
		-7.4033283286524815e-09,
		8.37169239652517e-11,
		-4.2796489422766765e-13
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx15(double t)
{
	static constexpr double P[] = {
		-2.0716173552085952,
		-0.6677865789134582,
		-0.9469183588702466,
		-0.003367793385122891,
		0.0001564020930025248,
		-5.19791945487372e-06,
		1.2039575142390604e-07,
		-1.8477994425430282e-09,
		1.6908573164121548e-11,
		-6.990670482225596e-14
	};

	double res = P[9];
	for (size_t i = 0; i < 9; i++)
		res = res * t + P[8 - i];
	return res;
}

static inline double Approx16(double t)
{
	static constexpr double P[] = {
		-2.8733800773456473,
		-0.42040398067011103,
		-0.9809615989177397,
		-0.0006207756095237024,
		1.2884643543110702e-05,
		-1.5262996617461687e-07,
		7.862698390437164e-10
	};

	double res = P[6];
	for (size_t i = 0; i < 6; i++)
		res = res * t + P[5 - i];
	return res;
}

static inline double AddEm(double x)
{
	static constexpr double emHigh = 0.36787944117144232160;
	static constexpr double emLow = -1.2428753672788363168e-17;

	return (x + emHigh) + emLow;
}

static inline double NearBranchWm1(double x)
{
	static constexpr double P[] = {
		-0.9999999999999998,
		-1.0000000000002183,
		-0.3333333332991474,
		-0.15277777988675093,
		-0.079629561585163,
		-0.04450363652960061,
		-0.025967945747804935,
		-0.015781789829297576,
		-0.008713980180658111,
		-0.010038019576077164,
		0.009199023794468941,
		-0.032793012029992144,
		0.04846364406273077,
		-0.056801153904189224,
		0.03763029838961123,
		-0.013204978244912841
	};

	static constexpr double s2e = 2.331643981597124;
	double p = sqrt(AddEm(x)) * s2e;

	// Evaluate polynomial using Horner's Method
	double value = P[15];
	for (size_t i = 0; i < 15; i++)
		value = value * p + P[14 - i];

	return value;
}

double MuirWm1v2(double x)
{
	if (x < -0.313486180883)
		return NearBranchWm1(x);

	double t = sqrt(-1 - log(-x));
	if (t < 6.23)
	{
		if (t < 2.35)
		{
			if (t < 1.22)
			{
				if (t < 0.79) return Approx1(t);
				return Approx2(t);
			}
			if (t < 1.72) return Approx3(t);
			return Approx4(t);
		}
		if (t < 4.05)
		{
			if (t < 3.18) return Approx5(t);
			return Approx6(t);
		}
		if (t < 5.05) return Approx7(t);
		return Approx8(t);
	}
	if (t < 13.98)
	{
		if (t < 9.34)
		{
			if (t < 7.64) return Approx9(t);
			return Approx10(t);
		}
		if (t < 11.42) return Approx11(t);
		return Approx12(t);
	}
	if (t < 21.22)
	{
		if (t < 17.18) return Approx13(t);
		return Approx14(t);
	}
	if (t < 26.37) return Approx15(t);
	return Approx16(t);
}
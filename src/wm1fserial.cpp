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
		0.00159381916122638080327842,-0.41846203646231866905097,13.0754476452340109046875,-33.425418804170640678635,-87.607345098466053529541
	};
	static constexpr double Q[] = {
		-0.000164817805247449945810037,0.062386189729457743417642,-2.90922590458749061427806,20.1957871204740373956932,1.
	};

	double t = (double)x;

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
		-12.09015688037941661084594,32135.56153372749530599325,-1.13911365007445329260617e7,6.600688192690400594562413e8,-2.136780678942179010830866e9
	};
	static constexpr double Q[] = {
		1.,-3474.789300680873162727,1.568926686034086152432803e6,-1.275251307979945845940914e8,1.261213777060697637244314e9
	};

	double t = (double)x;

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
		-14.56249821961777349366917,403320.3333527217959630796,-1.497879223894395170293063e9,9.62916113431681106029542e11,-6.03430015954642686243592e13
	};
	static constexpr double Q[] = {
		1.,-34265.87036329194996156926,1.521927455749773326092131e8,-1.22359941334666038804455e11,1.251300118362712875951786e13
	};

	double t = (double)x;

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
		-17.08650746361654600776381,5.222231087162755069389663e6,-2.131519875253351417813282e11,1.530697004513275159292723e15,-1.204303030667624689962588e18
	};
	static constexpr double Q[] = {
		1.,-364527.8285924629018789056,1.714851679428671852682649e10,-1.456699340067855240882402e14,1.583134020585771828059243e17
	};

	double t = (double)x;

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
		-19.65363217176511804577991,6.981171686805223179114112e7,-3.292362508083415294297803e13,2.747338886132200042460241e18,-2.632688582480525381859841e22
	};
	static constexpr double Q[] = {
		1.,-4.12717660940746078833663e6,2.189213015335675828455557e12,-2.089817431197110181778984e17,2.55141782750510018267205e21
	};

	double t = (double)x;

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
		-22.25842614929501931212822,9.624516363707524579641544e8,-5.492812670981184957840948e15,5.55616366190341123552896e21,-6.61117625096980985007073e26
	};
	static constexpr double Q[] = {
		1.,-4.926573687810429499529105e7,3.10794857350927823515454e14,-3.516710238245467777980717e20,5.079629764759798999192211e25
	};

	double t = (double)x;

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
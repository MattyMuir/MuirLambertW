#include "FukushimaMinimax.h"

#include <cmath>

#include "FukushimaMinimaxConstants.h"

template <typename Ty, size_t P, size_t Q>
static inline Ty EvaluateRational(Ty x, const Ty* coeff)
{
	Ty numer = coeff[P];
	for (size_t i = 0; i < P; i++)
		numer = numer * x + coeff[P - 1 - i];

	Ty denom = coeff[P + Q + 1];
	for (size_t i = 0; i < Q; i++)
		denom = denom * x + coeff[P + Q - i];

	return numer / denom;
}

template <size_t K>
static inline float RationalW0f(float x)
{
	double arg = (K < 18) ? sqrt(x + 0.36787944117144232160) : log((double)x);

	if constexpr (K == 1) return EvaluateRational<double, 4, 3>(arg, XUf);
	return EvaluateRational<double, 3, 3>(arg, XUf + 9 + (K - 2) * 8);
}

template <size_t K>
static inline double RationalW0(double x)
{
	double arg = (K < 18) ? sqrt(x + 0.36787944117144232160) : log(x);

	if constexpr (K == 1) return EvaluateRational<double, 8, 7>(arg, XU);
	return EvaluateRational<double, 7, 7>(arg, XU + 17 + (K - 2) * 16);
}

float FukushimaMinimaxW0(float x)
{
	if (x < Zk0f[0]) return RationalW0f<1>(x);
	if (x < Zk0f[1]) return RationalW0f<2>(x);
	if (x < Zk0f[2]) return RationalW0f<3>(x);
	if (x < Zk0f[3]) return RationalW0f<4>(x);
	if (x < Zk0f[4]) return RationalW0f<5>(x);
	if (x < Zk0f[5]) return RationalW0f<6>(x);
	if (x < Zk0f[6]) return RationalW0f<7>(x);
	if (x < Zk0f[7]) return RationalW0f<8>(x);
	if (x < Zk0f[8]) return RationalW0f<9>(x);
	if (x < Zk0f[9]) return RationalW0f<10>(x);
	if (x < Zk0f[10]) return RationalW0f<11>(x);
	if (x < Zk0f[11]) return RationalW0f<12>(x);
	if (x < Zk0f[12]) return RationalW0f<13>(x);
	if (x < Zk0f[13]) return RationalW0f<14>(x);
	if (x < Zk0f[14]) return RationalW0f<15>(x);
	if (x < Zk0f[15]) return RationalW0f<16>(x);
	if (x < Zk0f[16]) return RationalW0f<17>(x);
	return RationalW0f<18>(x);
}

double FukushimaMinimaxW0(double x)
{
	if (x < Zk0[0]) return RationalW0<1>(x);
	if (x < Zk0[1]) return RationalW0<2>(x);
	if (x < Zk0[2]) return RationalW0<3>(x);
	if (x < Zk0[3]) return RationalW0<4>(x);
	if (x < Zk0[4]) return RationalW0<5>(x);
	if (x < Zk0[5]) return RationalW0<6>(x);
	if (x < Zk0[6]) return RationalW0<7>(x);
	if (x < Zk0[7]) return RationalW0<8>(x);
	if (x < Zk0[8]) return RationalW0<9>(x);
	if (x < Zk0[9]) return RationalW0<10>(x);
	if (x < Zk0[10]) return RationalW0<11>(x);
	if (x < Zk0[11]) return RationalW0<12>(x);
	if (x < Zk0[12]) return RationalW0<13>(x);
	if (x < Zk0[13]) return RationalW0<14>(x);
	if (x < Zk0[14]) return RationalW0<15>(x);
	if (x < Zk0[15]) return RationalW0<16>(x);
	if (x < Zk0[16]) return RationalW0<17>(x);
	if (x < Zk0[17]) return RationalW0<18>(x);
	return RationalW0<19>(x);
}

float FukushimaMinimaxWm1(float x)
{
	return 0.0f;
}

double FukushimaMinimaxWm1(double x)
{
	return 0.0;
}
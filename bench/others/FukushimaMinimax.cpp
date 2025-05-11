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

float FukushimaMinimaxW0(float x)
{
	if (x < Zkf[0]) return RationalW0f<1>(x);
	if (x < Zkf[1]) return RationalW0f<2>(x);
	if (x < Zkf[2]) return RationalW0f<3>(x);
	if (x < Zkf[3]) return RationalW0f<4>(x);
	if (x < Zkf[4]) return RationalW0f<5>(x);
	if (x < Zkf[5]) return RationalW0f<6>(x);
	if (x < Zkf[6]) return RationalW0f<7>(x);
	if (x < Zkf[7]) return RationalW0f<8>(x);
	if (x < Zkf[8]) return RationalW0f<9>(x);
	if (x < Zkf[9]) return RationalW0f<10>(x);
	if (x < Zkf[10]) return RationalW0f<11>(x);
	if (x < Zkf[11]) return RationalW0f<12>(x);
	if (x < Zkf[12]) return RationalW0f<13>(x);
	if (x < Zkf[13]) return RationalW0f<14>(x);
	if (x < Zkf[14]) return RationalW0f<15>(x);
	if (x < Zkf[15]) return RationalW0f<16>(x);
	if (x < Zkf[16]) return RationalW0f<17>(x);
	return RationalW0f<18>(x);
}

double FukushimaMinimaxW0(double x)
{
	return 0.0;
}

float FukushimaMinimaxWm1(float x)
{
	return 0.0f;
}

double FukushimaMinimaxWm1(double x)
{
	return 0.0;
}
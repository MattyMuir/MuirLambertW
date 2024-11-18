#include "wm1fserial.h"

#include <cmath>
#include <cfenv>

// TODO - remove
#include <iostream>
#include <format>

float AddEm(float x)
{
	static constexpr float emHigh = 0.36787945f;
	static constexpr float emLow = -9.149756e-09f;
	return (x + emHigh) + emLow;
}

float IntPow(float f, uint64_t p)
{
	float ret = 1.0f;
	for (uint64_t i = 0; i < p; i++)
		ret *= f;
	return ret;
}

float MuirWm1(float x)
{
	float t;
	if (x < -0.2f)
		t = std::sqrt(std::log1p(AddEm(x) * -2.7182817f) * -2.0f);
	else
		t = std::sqrt(std::log(-x) * -2.0f - 2.0f);

	float P[] = {
		0.0f,
		-775.846669858888f,
		-869.712104557702f,
		-411.696629665943f,
		-100.601001722587f,
		-12.230662731677f,
		-0.4994929301872f,
		-4.04431805144266e-06f
	};

	float Q[] = {
		775.84618883946f,
		611.10071342493f,
		186.434267056385f,
		24.3695096706907f,
		1.0f
	};

	float numer = P[7];
	for (size_t i = 0; i < 7; i++)
		numer = numer * t + P[6 - i];

	float denom = Q[4];
	for (size_t i = 0; i < 4; i++)
		denom = denom * t + Q[3 - i];

	//std::cout << std::format("{}\n", numer);
	//std::cout << std::format("{}\n", denom);

	return numer / denom - 1.0f;
}
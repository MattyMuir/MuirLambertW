#include <iostream>
#include <format>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <MuirLambertW.h>

#define TIMER_NPRINT
#include "../bench/Timer.h"

static constexpr float EM_UPf = -0.36787942f;
static constexpr double EM_UP = -0.3678794411714423;

Intervalf ReferenceW0f(float x) { static ReferenceWf evaluator; return evaluator.W0(x); }
Intervalf ReferenceWm1f(float x) { static ReferenceWf evaluator; return evaluator.Wm1(x); }
Interval ReferenceW0(double x) { static ReferenceW evaluator; return evaluator.W0(x); }
Interval ReferenceWm1(double x) { static ReferenceW evaluator; return evaluator.Wm1(x); }

float ExpMapW0(float x)
{
	return EM_UPf + exp(x);
}

float ExpMapWm1(float x)
{
	if (x < 88.7f)
		return EM_UPf / (1 + expf(x));

	return EM_UPf / (1 + exp(x - 14)) * 8.315287e-07;
}

double ExpMapW0(double x)
{
	return EM_UP + exp(x);
}

double ExpMapWm1(double x)
{
	if (x < 700)
		return EM_UP / (1 + exp(x));

	return EM_UP / (1 + exp(x - 62)) * 1.185064864233981e-27;
}

int main()
{
	static std::mt19937_64 gen{ std::random_device{}() };
	static ReciprocalDistributionEx<float> dist{ EM_UPf, INFINITY, false };

	MaxULPRounded(ReferenceW0f, MakeSerial<float, MuirW0>, []() { return dist(gen); }, 0);
}
#include <iostream>
#include <charconv>
#include <random>

#include <ReferenceLambertW.h>
#include <flttestlib.h>
#include <MuirLambertW.h>

#define ERROR(msg) { std::cout << msg << '\n'; return 1; }

static constexpr float EM_UPf = -0.36787942f;
static constexpr double EM_UP = -0.3678794411714423;

Intervalf ReferenceW0f(float x) { static ReferenceWf evaluator; return evaluator.W0(x); }
Intervalf ReferenceWm1f(float x) { static ReferenceWf evaluator; return evaluator.Wm1(x); }
Interval ReferenceW0(double x) { static ReferenceW evaluator; return evaluator.W0(x); }
Interval ReferenceWm1(double x) { static ReferenceW evaluator; return evaluator.Wm1(x); }

template <typename Ty>
auto Reference(int64_t branch, Ty x)
{
	if constexpr (std::is_same_v<Ty, float>)
	{
		if (branch == 0) return ReferenceW0f(x);
		return ReferenceWm1f(x);
	}
	else
	{
		if (branch == 0) return ReferenceW0(x);
		return ReferenceWm1(x);
	}
}

template <typename Ty>
using FuncPtr = Ty(*)(Ty);

template <typename Ty, FuncPtr<Ty> Func>
Ty Overload(Ty x)
{
	return Func(x);
}

template <typename Ty> 
constexpr Ty GetEmUp()
{
	if (std::is_same_v<Ty, float>) return EM_UPf;
	return EM_UP;
}

template <typename Ty>
Ty GetMapMin(int64_t branch)
{
	if constexpr (std::is_same_v<Ty, float>)
	{
		return (branch == 0) ? -18.021828f : -16.635532f;
	}
	return (branch == 0) ? -38.123094930796995 : -36.73680056967711;
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

template <typename Ty>
int ExhaustiveTest(auto func, int64_t branch, UIntType<Ty> thresh)
{
	Ty min = GetEmUp<Ty>();
	Ty max = (branch == 0) ? INFINITY : 0;

	for (Ty x = min; x != max; x = std::nextafter(x, INFINITY))
	{
		Ty approx = func(x);
		auto exact = Reference(branch, x);
		if (ULPDistance(approx, exact) > thresh)
			ERROR(std::format("x: {}", x));
	}

	return 0;
}

template <typename Ty>
int UnitTests(auto func, int64_t branch)
{
	Ty emDown = std::nextafter(GetEmUp<Ty>(), -INFINITY);

	if (!std::isnan(func(emDown)))
		ERROR(std::format("x: {} W(x) should be NaN\n", emDown));

	if (!std::isnan(func(-INFINITY)))
		ERROR("W(-inf) should be NaN\n");

	static std::mt19937_64 gen{ std::random_device{}() };
	ReciprocalDistributionEx<Ty> beforeBranchDist{ -INFINITY, emDown, false };
	for (size_t i = 0; i < 100'000; i++)
	{
		Ty x = beforeBranchDist(gen);
		if (!std::isnan(func(x)))
			ERROR(std::format("x: {} W(x) should be NaN\n", x));
	}

	if (branch == 0)
	{
		if (func(0) != 0)
			ERROR("W(0) should equal 0");
		if (func(INFINITY) != INFINITY)
			ERROR("W(inf) should equal inf");
	}
	else
	{
		ReciprocalDistributionEx<Ty> posDist{ 0, INFINITY, false };
		for (size_t i = 0; i < 100'000; i++)
		{
			Ty x = posDist(gen);
			if (!std::isnan(func(x)))
				ERROR(std::format("x: {} W(x) should be NaN\n", x));
		}
	}

	return 0;
}

int DoRandomTests(auto func, int64_t branch, auto dist, size_t iter, uint64_t thresh, MapTy<double> map = IdentityMap)
{
	static std::mt19937_64 gen{ std::random_device{}() };
	
	for (size_t i = 0; i < iter; i++)
	{
		double x = map(dist(gen));
		double approx = func(x);
		auto exact = Reference(branch, x);

		if (ULPDistance(approx, exact) > thresh)
			ERROR(std::format("x: {}\n", x));
	}

	return 0;
}

int RandomTests(auto func, int64_t branch, size_t iter, uint64_t thresh)
{
	double min = GetEmUp<double>();
	double max = (branch == 0) ? INFINITY : 0;

	// Reciprocal test
	ReciprocalDistributionEx<double> dist1{ min, max, false };
	if (DoRandomTests(func, branch, dist1, iter, thresh)) return 1;

	// Near branch test
	std::uniform_real_distribution<double> dist2{ GetMapMin<double>(branch), -3.0 };
	if (DoRandomTests(func, branch, dist2, iter, thresh, (branch == 0) ? ExpMapW0 : ExpMapWm1)) return 1;

	// Problem area test
	std::uniform_real_distribution<double> dist3{ -0.35, -0.25 };
	if (DoRandomTests(func, branch, dist3, iter, thresh)) return 1;

	return 0;
}

int main(int argc, char** argv)
{
	// === Parameters ===
	static constexpr size_t Iter = 1'000'000;
	// ==================

#if 1
	// Check number of arguments is correct
	if (argc != 2)
		ERROR("Test must have exactly one extra argument");

	// Get second argument
	std::string arg{ argv[1] };

	// Convert to integer
	size_t testIdx;
	auto convRes = std::from_chars(arg.data(), arg.data() + arg.size(), testIdx);
	if (convRes.ec != std::errc())
		ERROR("Test index could not be parsed");
#else
	size_t testIdx = 0;
#endif

	switch (testIdx)
	{
	// Exhaustive tests
	case 0: return ExhaustiveTest<float>(Overload<float, MuirW0>, 0, 1);
	case 1: return ExhaustiveTest<float>(Overload<float, MuirWm1>, -1, 1);
	case 2: return ExhaustiveTest<float>(MakeSerial<float, MuirW0>, 0, 3);
	case 3: return ExhaustiveTest<float>(MakeSerial<float, MuirWm1>, -1, 4);

	// Unit tests
	case 4: return UnitTests<float>(Overload<float, MuirW0>, 0);
	case 5: return UnitTests<float>(Overload<float, MuirWm1>, -1);
	case 6: return UnitTests<float>(MakeSerial<float, MuirW0>, 0);
	case 7: return UnitTests<float>(MakeSerial<float, MuirWm1>, -1);
	case 8: return UnitTests<double>(Overload<double, MuirW0>, 0);
	case 9: return UnitTests<double>(Overload<double, MuirWm1>, -1);
	case 10: return UnitTests<double>(MakeSerial<double, MuirW0>, 0);
	case 11: return UnitTests<double>(MakeSerial<double, MuirWm1>, -1);

	// Sample tests
	case 12: return RandomTests(Overload<double, MuirW0>, 0, Iter, 4);
	case 13: return RandomTests(Overload<double, MuirWm1>, -1, Iter, 4);
	case 14: return RandomTests(MakeSerial<double, MuirW0>, 0, Iter, 4);
	case 15: return RandomTests(MakeSerial<double, MuirWm1>, -1, Iter, 4);
	}
}
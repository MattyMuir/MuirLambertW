#pragma once
#include "ReciprocalDistribution.h"

// Extended reciprocal distribution between any two endpoints
class ReciprocalDistributionEx
{
protected:
	enum class RangeType : uint8_t { JustNeg, PosNeg, Zero, JustPos };

public:
	ReciprocalDistributionEx() = default;
	ReciprocalDistributionEx(double min_, double max_, bool favorEndpoints_);

	template <typename Engine>
	double operator()(Engine& eng)
	{
		uint64_t branch;
		if (favorEndpoints)
		{
			branch = eng();
			if (branch < (uint64_t)((double)Engine::max() * 0.05))
				return min;
			if (branch > (uint64_t)((double)Engine::max() * 0.95))
				return max;
		}

		switch (rangeType)
		{
		case RangeType::JustNeg:
			return -negDist(eng);
		case RangeType::PosNeg:
			if (!favorEndpoints)
				branch = eng();
			if (branch % 2)
				return -negDist(eng);
			return posDist(eng);
		case RangeType::Zero:
			return 0.0;
		case RangeType::JustPos:
			return posDist(eng);
		}
	}

protected:
	double min, max;
	bool favorEndpoints;
	RangeType rangeType;
	ReciprocalDistribution negDist, posDist;

	static double SanitizeEndpoint(double val);
};
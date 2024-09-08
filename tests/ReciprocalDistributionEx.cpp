#include "ReciprocalDistributionEx.h"

static constexpr double DoubleMin = std::numeric_limits<double>::denorm_min();
static constexpr double DoubleMax = std::numeric_limits<double>::max();
static constexpr double Infinity = std::numeric_limits<double>::infinity();

ReciprocalDistributionEx::ReciprocalDistributionEx(double min_, double max_, bool favorEndpoints_)
	: min(min_), max(max_), favorEndpoints(favorEndpoints_)
{
	// Work out range type
	if (min < 0 && max > 0)
		rangeType = RangeType::PosNeg;
	else if (min < 0)
		rangeType = RangeType::JustNeg;
	else if (max > 0)
		rangeType = RangeType::JustPos;
	else
		rangeType = RangeType::Zero;

	// Sanitize endpoints
	double minSanitized = SanitizeEndpoint(min);
	double maxSanitized = SanitizeEndpoint(max);

	// Initialize child distributions
	switch (rangeType)
	{
	case RangeType::JustNeg:
		negDist = { maxSanitized, minSanitized };
		return;
	case RangeType::PosNeg:
	{
		negDist = { DoubleMin, minSanitized };
		posDist = { DoubleMin, maxSanitized };
		return;
	}
	case RangeType::JustPos:
		posDist = { minSanitized, maxSanitized };
		return;
	case RangeType::Zero:
		return;
	}
}

double ReciprocalDistributionEx::SanitizeEndpoint(double val)
{
	val = std::abs(val);

	if (val == 0)
		return DoubleMin;
	if (val == Infinity)
		return DoubleMax;

	return val;
}
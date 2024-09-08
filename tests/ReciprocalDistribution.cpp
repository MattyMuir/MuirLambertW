#include "ReciprocalDistribution.h"

#include <numeric>

ReciprocalDistribution::ReciprocalDistribution(double min_, double max_)
	: min(min_), max(max_), dist(0, 1)
{
	int minExp, maxExp;
	double minMan = std::frexp(min, &minExp);
	double maxMan = std::frexp(max, &maxExp);

	// Compute logWidth
	if (maxExp - minExp > 10)
		logWidth = std::log(max) - std::log(min);
	else
		logWidth = std::log1p((max - min) / min);

	// Compute geometricMean
	double meanMan = sqrt(minMan * maxMan);
	int expSum = minExp + maxExp;
	geometricMean = std::ldexp(meanMan, expSum / 2);
	if (expSum % 2)
		geometricMean *= 1.4142135623730950488;

	// Check if mean is needed
	double expMidpoint = exp(logWidth / 2);
	useMean = !std::isfinite(min * expMidpoint) || !std::isfinite(max / expMidpoint);
}

double ReciprocalDistribution::Min() const
{
	return min;
}

double ReciprocalDistribution::Max() const
{
	return max;
}
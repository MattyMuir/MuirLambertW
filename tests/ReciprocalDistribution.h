#pragma once
#include <random>

// Reciprocal distribution between two positive, finite endpoints
class ReciprocalDistribution
{
public:
	ReciprocalDistribution(double min_, double max_);

	double Min() const;
	double Max() const;

	template <typename Engine>
	double operator()(Engine& eng)
	{
		double val = dist(eng);

		if (useMean && val > 0.25 && val < 0.75)
			return geometricMean * std::exp((val - 0.5) * logWidth);

		if (val < 0.5)
			return min * std::exp(val * logWidth);
		return max / std::exp((1 - val) * logWidth);
	}

protected:
	double min, max, logWidth, geometricMean;
	bool useMean;
	std::uniform_real_distribution<double> dist;

	ReciprocalDistribution() = default;

	friend class ReciprocalDistributionEx;
};
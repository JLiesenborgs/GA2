#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"

namespace eatk
{

class RemainingTargetPopulationSizeIteration : public PopulationCrossoverIteration
{
public:
	RemainingTargetPopulationSizeIteration() : m_target(0) { }
	~RemainingTargetPopulationSizeIteration() { }
	
	void startNewIteration(const Population &newPopulation, size_t targetPopulationSize) override
	{
		m_target = targetPopulationSize;
	}
	
	bool iterate(const Population &newPopulation) override
	{
		if (newPopulation.size() >= m_target)
			return false;
		return true;
	}
private:
	size_t m_target;
};

}

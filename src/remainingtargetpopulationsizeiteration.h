#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"

namespace eatk
{

class RemainingTargetPopulationSizeIteration : public PopulationCrossoverIteration
{
public:
	RemainingTargetPopulationSizeIteration() : m_remaining(0) { }
	~RemainingTargetPopulationSizeIteration() { }
	
	void startNewIteration(const Population &newPopulation, size_t targetPopulationSize) override
	{
		if (targetPopulationSize < newPopulation.size())
			m_remaining = 0;
		else
			m_remaining = targetPopulationSize - newPopulation.size();
	}
	
	bool iterate(const Population &newPopulation) override
	{
		if (m_remaining == 0)
			return false;
		m_remaining--;
		return true;
	}
private:
	size_t m_remaining;
};

}

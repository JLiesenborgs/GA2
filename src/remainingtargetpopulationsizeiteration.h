#pragma once

#include "crossovermutation.h"

class RemainingTargetPopulationSizeIteration : public PopulationCrossoverIteration
{
public:
    RemainingTargetPopulationSizeIteration() : m_remaining(0) { }
    ~RemainingTargetPopulationSizeIteration() { }
	
    void startNewIteration(std::shared_ptr<Population> &newPopulation, size_t targetPopulationSize) override
    {
        if (targetPopulationSize < newPopulation->size())
            m_remaining = 0;
        else
            m_remaining = targetPopulationSize - newPopulation->size();
    }
	
	bool iterate(std::shared_ptr<Population> &newPopulation) override
    {
        if (m_remaining == 0)
            return false;
        m_remaining--;
        return true;
    }
private:
    size_t m_remaining;
};

#pragma once

#include "mogal2config.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace mogal2
{

template<class T>
class VectorGenomeUniformMutation : public GenomeMutation
{
public:
    VectorGenomeUniformMutation(double mutationFraction, T minValue, T maxValue, std::shared_ptr<RandomNumberGenerator> &rng)
        :  m_mutationFraction(mutationFraction), m_min(minValue), m_max(maxValue), m_rng(rng) { }
    ~VectorGenomeUniformMutation() { }

	errut::bool_t check(const Genome &genome) override
    {
        const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&genome);
        if (!pGenome)
            return "Genome is not of the expected type";
        return true;
    }

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
    {
        VectorGenome<T> &g = static_cast<VectorGenome<T>&>(genome);
        std::vector<T> &v = g.getValues();
    
        for (auto &x : v)
        {
            // TODO: can we do this without as many random numbers? Floats? Generate indices first by some other means?
            if (m_rng->getRandomDouble() < m_mutationFraction)
            {
                x = m_rng->getRandomDouble(m_min, m_max);
                isChanged = true;
            }
        }
        return true;
    }
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    double m_mutationFraction;
    double m_min, m_max;
};

}

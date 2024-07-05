#pragma once

#include "eatkconfig.h"
#include "vectorgenomefitness.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include <cassert>

namespace eatk
{

template<class T, class GS>
class SimulatedBinaryCrossoverTemplate : public GenomeCrossover
{
public:
    SimulatedBinaryCrossoverTemplate(const std::shared_ptr<RandomNumberGenerator> &rng, T n) 
        : GenomeCrossover(2), m_rng(rng), m_n(n) // Two parents
    {
        m_exponent = (T)(1.0/(m_n+1.0));
    }

	errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
    {
        if (parents.size() != 2)
            return "Expecting two parents";

        for (auto &g : parents)
		{
			const GS *pGenome = dynamic_cast<const GS *>(g.get());
			if (!pGenome) 
				return "Parents contain an unexpected type: " + std::string(typeid(*pGenome).name());
		}

        return true;
    }

	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
                                    std::vector<std::shared_ptr<Genome>> &generatedOffspring) override
    {
        assert(parents.size() == 2);
        
        T p = (T)m_rng->getRandomFloat();
        bool lessThanHalf = true;
        if (p > (T)0.5)
        {
            lessThanHalf = false;
            p = ((T)1)-p;
        }
        T X = std::pow(((T)2)*p, m_exponent);
        T beta = (lessThanHalf)?X:((T)1)/X;
                                  
        std::shared_ptr<Genome> c[2] = { parents[0]->createCopy(false), parents[0]->createCopy(false) };
        const Genome &p1 = *parents[0];
        const Genome &p2 = *parents[1];

        auto createChild = [&p1, &p2](Genome &c, T beta)
        {
            size_t N = VectorGenome<T>::getSize(c);
            for (size_t i = 0 ; i < N ; i++)
            {
                T val = ((T)0.5) * ( VectorGenome<T>::getValue(p1, i)*(((T)1)+beta)
                                   + VectorGenome<T>::getValue(p2, i)*(((T)1)-beta) );
                VectorGenome<T>::setValue(c, i, val);
            }
        };

        createChild(*c[0], beta);
        createChild(*c[1], -beta);

        generatedOffspring.clear();
        generatedOffspring.push_back(c[0]);
        generatedOffspring.push_back(c[1]);

        return true;
    }
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    T m_n, m_exponent;
};

template <class T>
using VectorGenomeSimulatedBinaryCrossover = SimulatedBinaryCrossoverTemplate<T, VectorGenome<T>>;


}
#pragma once

#include "eatkconfig.h"
#include "vectorgenomefitness.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include <cassert>

namespace eatk
{

template<class T>
class VectorGenomeSimulatedBinaryCrossover : public GenomeCrossover
{
public:
    VectorGenomeSimulatedBinaryCrossover(const std::shared_ptr<RandomNumberGenerator> &rng, 
                                         T n) 
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
			const Genome *pGenome = g.get();
			const VectorGenome<T> *pVectorGenome = dynamic_cast<const VectorGenome<T> *>(pGenome);
			if (!pVectorGenome) 
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
            p = (T)1-p;
        }
        T X = std::pow(((T)2)*p, m_exponent);
        T beta = (lessThanHalf)?X:((T)1)/X;

        assert(dynamic_cast<VectorGenome<T>*>(parents[0].get()) && dynamic_cast<VectorGenome<T>*>(parents[1].get()));
        VectorGenome<T> &p1 = static_cast<VectorGenome<T>&>(*parents[0]);
        VectorGenome<T> &p2 = static_cast<VectorGenome<T>&>(*parents[1]);
                                  
        std::shared_ptr<Genome> c[2] = { p1.createCopy(false), p1.createCopy(false) };
        
        auto createChild = [&p1, &p2](Genome &c, T beta)
        {
            for (size_t i = 0 ; i < VectorGenome<T>::getSize(c) ; i++)
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

}
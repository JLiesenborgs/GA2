#pragma once

#include "mogal2config.h"
#include "vectorgenomefitness.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace mogal2
{

template<class T>
class UniformVectorGenomeCrossover : public GenomeCrossover
{
public:
    UniformVectorGenomeCrossover(const std::shared_ptr<RandomNumberGenerator> &rng, bool twoOffspring = false) 
        : m_rng(rng), m_twoOffspring(twoOffspring) { }
    ~UniformVectorGenomeCrossover() { }
	
    errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
    {
        if (parents.size() != 2)
            return "Expecting two parents";

        int numAlleles = -1;
        for (auto &g : parents)
        {
            const Genome *pGenome = g.get();
            const VectorGenome<T> *pVectorGenome = dynamic_cast<const VectorGenome<T> *>(pGenome);
            if (!pVectorGenome) 
                return "Parents contain an unexpected type: " + std::string(typeid(*pGenome).name());

            int n = (int)pVectorGenome->getValues().size();
            if (numAlleles >= 0 && n != numAlleles)
                return "Different number of alleles in parents";
            numAlleles = n;
        }
        return true;
    }

	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
                                    std::vector<std::shared_ptr<Genome>> &offspring) override
    {
        // TODO: just an assert instead?
        if (parents.size() != 2)
            return "Expecting two parents";

        const VectorGenome<T> *pParents[2];
        for (int i = 0 ; i < 2 ; i++)
        {
            pParents[i] = static_cast<const VectorGenome<T>*>(parents[i].get());
            assert(pParents[i]);
        }

        const std::vector<T> *values[2] = { &pParents[0]->getValues(), &pParents[1]->getValues() };
        
        const int numOffspring = (m_twoOffspring)?2:1;
        offspring.resize(numOffspring);

        // TODO: create something that recycles the genomes? Need to watch out for overwriting existing ones though!!
        for (auto &o : offspring)
            o = pParents[0]->createCopy(false);

        // TODO: instead of using floats per allele, generate 32bit numbers and use the bits
        //       to decide which parent
        auto pickParent = [this]()
        {
            float x = m_rng->getRandomFloat();
            return (x < 0.5f)?0:1;
        };

        size_t numAlleles = values[0]->size();
        for (size_t i = 0 ; i < numAlleles ; i++)
        {
            int p1 = pickParent();
            int p2 = 1-p1;

            auto setAllele = [i,values,offspring](int offspringIndex, int parent)
            {
                std::vector<T> &o = static_cast<VectorGenome<T> &>(*offspring[offspringIndex]).getValues();
                o[i] = (*values[parent])[i];    
            };
            setAllele(0, p1);
            if(m_twoOffspring)
                setAllele(1, p2);
        }
        return true;
    }
private:
    std::shared_ptr<RandomNumberGenerator> m_rng;
    const bool m_twoOffspring;
};

}

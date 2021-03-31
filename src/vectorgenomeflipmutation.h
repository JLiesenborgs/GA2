#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "vectorgenomefitness.h"

namespace eatk
{

template<class T>
class VectorGenomeFlipMutation : public GenomeMutation
{
public:
	VectorGenomeFlipMutation(const std::shared_ptr<RandomNumberGenerator> &rng, double prob)
	 : m_rng(rng), m_prob(prob)
	{
	}

	errut::bool_t check(const Genome &genome) override
	{
		auto g = dynamic_cast<const VectorGenome<T> *>(&genome);
		if (!g)
			return "Genome type is incompatible";
		return true;
	}

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
	{
		VectorGenome<T> &g = static_cast<VectorGenome<T> &>(genome);
		auto &v = g.getValues();

		for (auto &x : v)
		{
			if (m_rng->getRandomDouble() < m_prob)
				x = !x;
		}
		return true;
	}
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	double m_prob;
};

}
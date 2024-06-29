#pragma once

#include "eatkconfig.h"
#include "vectorgenomefitness.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include <memory>
#include <limits>
#include <cmath>

namespace eatk
{

template <class T>
class VectorGenomeDELikeCrossOver : public GenomeCrossover
{
public:
	VectorGenomeDELikeCrossOver(std::shared_ptr<RandomNumberGenerator> &rng,
	                            bool extraParent = true,
	                            float F = std::numeric_limits<float>::quiet_NaN(), // NaN signals uniform number between 0 and 1 each time
	                            float CR = std::numeric_limits<float>::quiet_NaN() // same
	                            )
		: GenomeCrossover((extraParent)?4:3),
		  m_rng(rng),
		  m_F(F),
		  m_CR(CR),
		  m_extraParent(extraParent)
	{
	}

	errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
	{
		for (auto &p : parents)
		{
			if (!dynamic_cast<VectorGenome<T>*>(p.get()))
				return "Parent is of wrong type";
		}
		return true;
	}

	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
									std::vector<std::shared_ptr<Genome>> &generatedOffspring) override
	{
		assert((m_extraParent && parents.size() == 4) || (!m_extraParent && parents.size() == 3));
		generatedOffspring.clear();

		std::vector<T> *vg[3];
		for (size_t i = 0 ; i < 3 ; i++)
		{
			assert(dynamic_cast<VectorGenome<T>*>(parents[i].get()));
			VectorGenome<T> *pG = static_cast<VectorGenome<T>*>(parents[i].get());
			vg[i] = &pG->getValues();
		}

		auto result = parents[0]->createCopy();
		std::vector<T> &vo = (static_cast<VectorGenome<T>&>(*result)).getValues();

		size_t refParent = 3;
		if (!m_extraParent)
		{
			// Using 0 as reference parent doesn't seem to work all that well
			// Then again, that's also the one we're trying to adjust based on
			// the difference of the other two
			if (m_rng->getRandomFloat() < 0.5f)
				refParent = 1;
			else
				refParent = 2;
		}

		float CR = (isnan(m_CR))?m_rng->getRandomFloat():m_CR;
		float F = (isnan(m_F))?m_rng->getRandomFloat():m_F;

		std::vector<T> &c = static_cast<VectorGenome<T>&>(*parents[refParent]).getValues();
		size_t rndIdx = (size_t)m_rng->getRandomUint32() % vo.size();

		for (size_t i = 0 ; i < vo.size() ; i++)
		{
			if (m_rng->getRandomFloat() < CR || i == rndIdx)
				vo[i] = (*(vg[0]))[i] + F * ((*(vg[2]))[i] - (*(vg[1]))[i]);
			else
				vo[i] = c[i];
		}

		generatedOffspring.push_back(result);
		return true;
	}
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	float m_F, m_CR;
	bool m_extraParent;
};

}
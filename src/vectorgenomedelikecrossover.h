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

template <class T, class GS> // T is for data type (float/double), GS is to get/set vector value
class DELikeCrossOverTemplate : public GenomeCrossover
{
public:
	DELikeCrossOverTemplate(std::shared_ptr<RandomNumberGenerator> &rng,
	                            bool extraParent = true,
	                            float F = std::numeric_limits<float>::quiet_NaN(), // NaN signals uniform number between 0 and 1 each time
	                            float CR = std::numeric_limits<float>::quiet_NaN(), // same
	                            const std::vector<T> &lowerBound = std::vector<T>(),
								const std::vector<T> &upperBound = std::vector<T>())
		: GenomeCrossover((extraParent)?4:3),
		  m_rng(rng),
		  m_F(F),
		  m_CR(CR),
		  m_extraParent(extraParent),
		  m_lowerBound(lowerBound),
		  m_upperBound(upperBound)
	{
	}

	errut::bool_t check(const std::vector<std::shared_ptr<Genome>> &parents) override
	{
		if (parents.size() < 1)
			return "At least one parent must be present in check";

		for (auto &p : parents)
		{
			if (!dynamic_cast<VectorGenome<T>*>(p.get()))
				return "Parent is of wrong type";
		}

		if (m_upperBound.size() != 0)
		{
			if (m_upperBound.size() != GS::getSize(*parents[0]))
				return "Incorrect upperbound size";
		}

		if (m_lowerBound.size() != 0)
		{
			if (m_lowerBound.size() != GS::getSize(*parents[0]))
				return "Incorrect lowerbound size";
		}

		if (m_lowerBound.size() > 0 && m_upperBound.size() > 0)
		{
			for (size_t i = 0 ; i < m_lowerBound.size() ; i++)
				if (m_lowerBound[i] >= m_upperBound[i])
					return "Lower bound must be strictly less than upper bound everywhere";
		}

		return true;
	}

	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<Genome>> &parents,
									std::vector<std::shared_ptr<Genome>> &generatedOffspring) override
	{
		assert((m_extraParent && parents.size() == 4) || (!m_extraParent && parents.size() == 3));
		generatedOffspring.clear();

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

		auto result = parents[0]->createCopy();
		float CR = (isnan(m_CR))?m_rng->getRandomFloat():m_CR;
		float F = (isnan(m_F))?m_rng->getRandomFloat():m_F;

		size_t dimensions = GS::getSize(*result);
		size_t rndIdx = (size_t)m_rng->getRandomUint32() % dimensions;
		bool useLowerBound = (m_lowerBound.size() > 0);
		bool useUpperBound = (m_upperBound.size() > 0);

		for (size_t i = 0 ; i < dimensions ; i++)
		{
			T val = 0;
			if (m_rng->getRandomFloat() < CR || i == rndIdx)
				val = GS::getValue(*parents[0], i) + F * (GS::getValue(*parents[2], i) - GS::getValue(*parents[1], i));
			else
				val = GS::getValue(*parents[refParent],i);

			if (useLowerBound && val < m_lowerBound[i])
				val = ( m_lowerBound[i] + GS::getValue(*parents[refParent],i) )/((T)2.0);
			if (useUpperBound && val > m_upperBound[i])
				val = ( m_upperBound[i] + GS::getValue(*parents[refParent],i) )/((T)2.0);

			GS::setValue(*result, i, val);
		}

		generatedOffspring.push_back(result);
		return true;
	}
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	float m_F, m_CR;
	bool m_extraParent;
	std::vector<T> m_lowerBound, m_upperBound;
};

template <class T>
using VectorGenomeDELikeCrossOver = DELikeCrossOverTemplate<T, VectorGenome<T>>;

}
#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace eatk
{

template<class T>
class VectorGenomeUniformMutation : public GenomeMutation
{
public:
	VectorGenomeUniformMutation(double mutationFraction, T minValue, T maxValue, const std::shared_ptr<RandomNumberGenerator> &rng)
		:  m_mutationFraction(mutationFraction),
		   m_minValues(std::vector<T>(1, minValue)),
		   m_maxValues(std::vector<T>(1, maxValue)), m_rng(rng) { }

	VectorGenomeUniformMutation(double mutationFraction,
	                            const std::vector<T> &minValues,
								const std::vector<T> &maxValues,
								const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_mutationFraction(mutationFraction), m_minValues(minValues), m_maxValues(maxValues), m_rng(rng) { }
	~VectorGenomeUniformMutation() { }

	errut::bool_t check(const Genome &genome) override
	{
		const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&genome);
		if (!pGenome)
			return "Genome is not of the expected type";

		size_t s = pGenome->getValues().size();
		if (s == 0)
			return "Genome has no values";

		auto resizeToGenomeLength = [s](std::vector<T> &v) -> errut::bool_t
		{
			if (v.size() == 0)
				return "No min or max values set";

			if (v.size() > 1)
			{
				if (v.size() != s)
					return "Genome length is incompatible with min/max length";
				return true;
			}

			T value = v[0];
			v.resize(s, value);
			return true;
		};

		errut::bool_t r;

		if (!(r = resizeToGenomeLength(m_minValues)) || !(r = resizeToGenomeLength(m_maxValues)))
			return r;

		return true;
	}

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
	{
		VectorGenome<T> &g = static_cast<VectorGenome<T>&>(genome);
		std::vector<T> &v = g.getValues();

		assert(m_minValues.size() == v.size() && m_maxValues.size() == v.size());
	
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			// TODO: can we do this without as many random numbers? Floats? Generate indices first by some other means?
			if (m_rng->getRandomDouble() < m_mutationFraction)
			{
				v[i] = m_rng->getRandomDouble(m_minValues[i], m_maxValues[i]);
				isChanged = true;
			}
		}
		return true;
	}
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	double m_mutationFraction;
	std::vector<T> m_minValues, m_maxValues;
};

}

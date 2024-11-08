#pragma once

#include "eatkconfig.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"

namespace eatk
{

template<class T>
class VectorGenomeMutationBase : public GenomeMutation
{
public:
	VectorGenomeMutationBase(double mutationFraction, T minValue, T maxValue, const std::shared_ptr<RandomNumberGenerator> &rng)
		:  m_mutationFraction(mutationFraction),
		   m_minValues(std::vector<T>(1, minValue)),
		   m_maxValues(std::vector<T>(1, maxValue)), m_rng(rng) { }

	VectorGenomeMutationBase(double mutationFraction,
	                            const std::vector<T> &minValues,
								const std::vector<T> &maxValues,
								const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_mutationFraction(mutationFraction), m_minValues(minValues), m_maxValues(maxValues), m_rng(rng) { }
	~VectorGenomeMutationBase() { }

	errut::bool_t check(const Genome &genome) override
	{
		const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&genome);
		if (!pGenome)
			return "Genome is not of the expected type";

		size_t s = pGenome->getValues().size();
		if (s == 0)
			return "Genome has no values";

		for (size_t i = 0 ; i < m_minValues.size() ; i++)
			if (m_minValues[i] >= m_maxValues[i])
				return "Minimum value is not smaller than maximum at position " + std::to_string(i);

		errut::bool_t r;

		if (!(r = resizeToGenomeLength(s, m_minValues, "minimum")) || !(r = resizeToGenomeLength(s, m_maxValues, "maximum")))
			return r;

		return true;
	}

protected:
	static errut::bool_t resizeToGenomeLength(size_t s, std::vector<T> &v, const std::string &arrayTypeName)
	{
		if (v.size() == 0)
			return "No " + arrayTypeName + " values set";

		if (v.size() > 1)
		{
			if (v.size() != s)
				return "Genome length is incompatible with " + arrayTypeName + " values array length";
			return true;
		}

		T value = v[0];
		v.resize(s, value);
		return true;
	}

	std::shared_ptr<RandomNumberGenerator> m_rng;
	double m_mutationFraction;
	std::vector<T> m_minValues, m_maxValues;
};

template<class T>
class VectorGenomeUniformMutation : public VectorGenomeMutationBase<T>
{
public:
	VectorGenomeUniformMutation(double mutationFraction, T minValue, T maxValue, const std::shared_ptr<RandomNumberGenerator> &rng)
		: VectorGenomeMutationBase<T>(mutationFraction, minValue, maxValue, rng) { }

	VectorGenomeUniformMutation(double mutationFraction,
	                            const std::vector<T> &minValues,
								const std::vector<T> &maxValues,
								const std::shared_ptr<RandomNumberGenerator> &rng)
		: VectorGenomeMutationBase<T>(mutationFraction, minValues, maxValues, rng) { }

	~VectorGenomeUniformMutation() { }

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
	{
		VectorGenome<T> &g = static_cast<VectorGenome<T>&>(genome);
		std::vector<T> &v = g.getValues();

		assert(VectorGenomeMutationBase<T>::m_minValues.size() == v.size() &&
		       VectorGenomeMutationBase<T>::m_maxValues.size() == v.size());
	
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			// TODO: can we do this without as many random numbers? Floats? Generate indices first by some other means?
			if (VectorGenomeMutationBase<T>::m_rng->getRandomDouble() < VectorGenomeMutationBase<T>::m_mutationFraction)
			{
				v[i] = VectorGenomeMutationBase<T>::m_rng->getRandomDouble(VectorGenomeMutationBase<T>::m_minValues[i],
																		   VectorGenomeMutationBase<T>::m_maxValues[i]);
				isChanged = true;
			}
		}
		return true;
	}
};

}

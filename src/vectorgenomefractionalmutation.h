#pragma once

#include "eatkconfig.h"
#include "vectorgenomeuniformmutation.h" // For base class
#include <limits>

namespace eatk
{

template<class T, class FractionalAdjuster, bool useReferenceSizes>
class VectorGenomeFractionalMutation : public VectorGenomeMutationBase<T>
{
public:
	VectorGenomeFractionalMutation(double mutationFraction, const std::shared_ptr<RandomNumberGenerator> &rng,
								   T hardMinValue = -std::numeric_limits<T>::infinity(),
								   T hardMaxValue = std::numeric_limits<T>::infinity())
		: VectorGenomeMutationBase<T>(mutationFraction, hardMinValue, hardMaxValue, rng) { }

	VectorGenomeFractionalMutation(double mutationFraction, T referenceSize,
								   const std::shared_ptr<RandomNumberGenerator> &rng,
								   T hardMinValue = -std::numeric_limits<T>::infinity(),
								   T hardMaxValue = std::numeric_limits<T>::infinity())
		: VectorGenomeMutationBase<T>(mutationFraction, hardMinValue, hardMaxValue, rng),
	      m_referenceSizes(std::vector<T>(1, referenceSize)) { }

	VectorGenomeFractionalMutation(double mutationFraction,
	                            const std::vector<T> &hardMinValues,
								const std::vector<T> &hardMaxValues,
								const std::shared_ptr<RandomNumberGenerator> &rng)
		: VectorGenomeMutationBase<T>(mutationFraction, hardMinValues, hardMaxValues, rng) { }

	VectorGenomeFractionalMutation(double mutationFraction,
								const std::vector<T> &referenceSizes,
	                            const std::vector<T> &hardMinValues,
								const std::vector<T> &hardMaxValues,
								const std::shared_ptr<RandomNumberGenerator> &rng)
		: VectorGenomeMutationBase<T>(mutationFraction, hardMinValues, hardMaxValues, rng),
	      m_referenceSizes(referenceSizes) { }

	~VectorGenomeFractionalMutation() { }

	errut::bool_t check(const Genome &genome) override
	{
		errut::bool_t r = VectorGenomeMutationBase<T>::check(genome);
		if (!r)
			return r;

		const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&genome);
		if (!pGenome)
			return "Genome is not of the expected type";

		size_t s = pGenome->getValues().size();
		if (useReferenceSizes)
		{
			if (!(r = resizeToGenomeLength(s, m_referenceSizes)))
				return r;
		}
		else
		{
			if (m_referenceSizes.size() != 0)
				return "Reference sizes set but not using them";
		}

		return true;
	}

	errut::bool_t mutate(Genome &genome, bool &isChanged) override
	{
		const std::vector<T> &hardMinValues = VectorGenomeMutationBase<T>::m_minValues;
		const std::vector<T> &hardMaxValues = VectorGenomeMutationBase<T>::m_maxValues;

		VectorGenome<T> &g = static_cast<VectorGenome<T>&>(genome);
		std::vector<T> &v = g.getValues();

		assert(hardMinValues.size() == v.size() && hardMaxValues.size() == v.size());
	
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			// TODO: can we do this without as many random numbers? Floats? Generate indices first by some other means?
			if (VectorGenomeMutationBase<T>::m_rng->getRandomDouble() < VectorGenomeMutationBase<T>::m_mutationFraction)
			{
				T oldValue = v[i];
				assert(oldValue >= hardMinValues[i] && oldValue <= hardMaxValues[i]);

				T newValue;
				if (!useReferenceSizes)
					newValue = FractionalAdjuster::adjustValue(oldValue, hardMinValues[i], hardMaxValues[i], *VectorGenomeMutationBase<T>::m_rng);
				else
					newValue = FractionalAdjuster::adjustValue(oldValue, m_referenceSizes[i], hardMinValues[i], hardMaxValues[i], *VectorGenomeMutationBase<T>::m_rng);

				if (newValue < hardMinValues[i])
					newValue = hardMinValues[i] + (oldValue - hardMinValues[i])/(T)2;

				if (newValue > hardMaxValues[i])
					newValue = hardMaxValues[i] - (hardMaxValues[i] - oldValue)/(T)2;

				v[i] = newValue;
				isChanged = (oldValue != newValue);
			}
		}
		return true;
	}

private:
	std::vector<T> m_referenceSizes;
};

}

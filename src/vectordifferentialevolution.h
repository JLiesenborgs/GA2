#pragma once

#include "eatkconfig.h"
#include "differentialevolutionevolver.h"
#include "vectorgenomefitness.h"
#include "valuefitness.h"



#include <iostream>

namespace eatk
{

template<class T>
class VectorDifferentialEvolutionMutation : public DifferentialEvolutionMutation
{
public:
	errut::bool_t check(const Genome &g) override;
	std::shared_ptr<Genome> mutate(const std::vector<const Genome*> &genomes, const std::vector<double> &weights) override;
};

template<class T>
class VectorDifferentialEvolutionCrossover : public DifferentialEvolutionCrossover
{
public:
	VectorDifferentialEvolutionCrossover(const std::shared_ptr<RandomNumberGenerator> &rng,
										 const std::vector<T> &lowerBounds = {}, const std::vector<T> &upperBounds = {})
		: m_rng(rng), m_lower(lowerBounds), m_upper(upperBounds)
	{
		m_useBounds = (m_lower.size()||m_upper.size() > 0)?true:false;
	}

	errut::bool_t check(const Genome &g) override;
	errut::bool_t crossover(double CR, Genome &mutantDest, const Genome &origVector) override;
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::vector<T> m_lower, m_upper;
	bool m_useBounds;
};

template<class VT, class FT>
class VectorDifferentialEvolutionIndividualCreation : public IndividualCreation
{
public:
	VectorDifferentialEvolutionIndividualCreation(const std::vector<VT> &lower, const std::vector<VT> &upper,
	                     						  const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_lower(lower), m_upper(upper), m_rng(rng) { 	}

	std::shared_ptr<Genome> createInitializedGenome() override;
	std::shared_ptr<Genome> createUnInitializedGenome() override;
	std::shared_ptr<Fitness> createEmptyFitness() override { return std::make_shared<ValueFitness<FT>>(); }
private:
	std::vector<VT> m_lower, m_upper;
	std::shared_ptr<RandomNumberGenerator> m_rng;
};

template<class T>
inline errut::bool_t VectorDifferentialEvolutionMutation<T>::check(const Genome &g)
{
	const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&g);
	if (!pGenome)
		return "Genome is not a VectorGenome of the correct type";
	return true;
}

template<class T>
inline std::shared_ptr<Genome> VectorDifferentialEvolutionMutation<T>::mutate(const std::vector<const Genome*> &genomes, const std::vector<double> &weights)
{
	assert(genomes.size() > 0);
	assert(genomes.size() == weights.size());

	std::shared_ptr<Genome> result = genomes[0]->createCopy();
	VectorGenome<T> &g0 = static_cast<VectorGenome<T>&>(*result);
	std::vector<T> &v0 = g0.getValues();

	v0.assign(v0.size(), 0);

	for (size_t j = 0 ; j < genomes.size() ; j++)
	{
		T f = (T)weights[j];

		const VectorGenome<T> &g1 = static_cast<const VectorGenome<T>&>(*genomes[j]);
		const std::vector<T> &v1 = g1.getValues();
	
		assert(v1.size() == v0.size());
		
		for (size_t i = 0 ; i < v0.size() ; i++)
			v0[i] += f * v1[i];
	}

	return result;
}

template<class T>
inline errut::bool_t VectorDifferentialEvolutionCrossover<T>::check(const Genome &g)
{
	const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&g);
	if (!pGenome)
		return "Genome is not a VectorGenome of the correct type";

	if (m_useBounds)
	{
		if (m_upper.size() != m_lower.size())
			return "Specified bound have different sizes";
		for (size_t i = 0 ; i < m_lower.size() ; i++)
		{
			if (m_lower[i] >= m_upper[i])
				return "Lower bound " + std::to_string(i) + " is larger or equal than upper bound";
		}
	}
	return true;
}

template<class T>
inline errut::bool_t VectorDifferentialEvolutionCrossover<T>::crossover(double CR, Genome &mutantDest, const Genome &origVector)
{
	VectorGenome<T> &g0 = static_cast<VectorGenome<T>&>(mutantDest);
	const VectorGenome<T> &g1 = static_cast<const VectorGenome<T>&>(origVector);

	std::vector<T> &v0 = g0.getValues();
	const std::vector<T> &v1 = g1.getValues();
	assert(v0.size() == v1.size());

	size_t rndIdx = ((size_t)m_rng->getRandomUint32())%(v0.size());
	for (size_t i = 0 ; i < v0.size() ; i++)
	{
		if (i != rndIdx && m_rng->getRandomDouble() > CR)
			v0[i] = v1[i];
	}

	if (m_useBounds) // v1 should already be within bounds, we won't check this (only in assert)!
	{
		for (size_t i = 0 ; i < v0.size() ; i++)
		{
			if (v0[i] < m_lower[i])
				v0[i] = (m_lower[i] + v1[i])/2.0;
			
			if (v0[i] > m_upper[i])
				v0[i] = (m_upper[i] + v1[i])/2.0;

			assert(v0[i] >= m_lower[i] && v0[i] <= m_upper[i]);
			assert(v1[i] >= m_lower[i] && v1[i] <= m_upper[i]);
		}
	}
	return true;
}

template<class VT, class FT>
inline std::shared_ptr<Genome> VectorDifferentialEvolutionIndividualCreation<VT,FT>::createInitializedGenome()
{
	if (m_lower.size() != m_upper.size())
		return nullptr;
	if (m_lower.size() == 0)
		return nullptr;
	
	std::shared_ptr<VectorGenome<VT>> genome = std::make_shared<VectorGenome<VT>>(m_lower.size());
	std::vector<VT> &values = genome->getValues();
	for (size_t i = 0 ; i < values.size() ; i++)
	{
		VT r = (VT)m_rng->getRandomDouble();
		values[i] = m_lower[i] + r*(m_upper[i] - m_lower[i]);
	}
	
	return genome;
}

template<class VT, class FT>
inline std::shared_ptr<Genome> VectorDifferentialEvolutionIndividualCreation<VT,FT>::createUnInitializedGenome()
{
	if (m_lower.size() != m_upper.size())
		return nullptr;
	if (m_lower.size() == 0)
		return nullptr;
	
	std::shared_ptr<VectorGenome<VT>> genome = std::make_shared<VectorGenome<VT>>(m_lower.size());
	std::vector<VT> &values = genome->getValues();
	for (size_t i = 0 ; i < values.size() ; i++)
		values[i] = std::numeric_limits<VT>::quiet_NaN();
	
	return genome;
}

}

#pragma once

#include "eatkconfig.h"
#include "differentialevolutionevolver.h"
#include "vectorgenomefitness.h"
#include "valuefitness.h"

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
	VectorDifferentialEvolutionCrossover(const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_rng(rng) { }

	errut::bool_t check(const Genome &g) override;
	errut::bool_t crossover(double CR, Genome &mutantDest, const Genome &origVector) override;
private:
	std::shared_ptr<RandomNumberGenerator> m_rng;
};

template<class VT, class FT>
class VectorDifferentialEvolutionIndividualCreation : public IndividualCreation
{
public:
	VectorDifferentialEvolutionIndividualCreation(const std::vector<VT> &lower, const std::vector<VT> &upper,
	                     						  const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_lower(lower), m_upper(upper), m_rng(rng) { 	}

	std::shared_ptr<Genome> createInitializedGenome() override;
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

#if 1
	size_t rndIdx = ((size_t)m_rng->getRandomUint32())%(v0.size());
	for (size_t i = 0 ; i < v0.size() ; i++)
	{
		if (i != rndIdx && m_rng->getRandomDouble() > CR)
			v0[i] = v1[i];
	}
#else
	size_t j = (size_t)m_rng->getRandomUint32()%(v0.size());
	for (size_t k = 1 ; k <= v0.size() ; k++)
	{
		if (k != v0.size() && m_rng->getRandomDouble() > CR)
			v0[j] = v1[j];

		j = (j+1)%v0.size();
	}
#endif
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


}

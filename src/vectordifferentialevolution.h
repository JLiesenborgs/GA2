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
	VectorDifferentialEvolutionMutation(T factor) : m_factor(factor) { }

	errut::bool_t check(const Genome &g) override;
	std::shared_ptr<Genome> mutate(const Genome &r1, const Genome &r2, const Genome &r3) override;
private:
	const T m_factor;
};

template<class T>
class VectorDifferentialEvolutionCrossover : public DifferentialEvolutionCrossover
{
public:
	VectorDifferentialEvolutionCrossover(double CR, const std::shared_ptr<RandomNumberGenerator> &rng)
		: m_CR(CR), m_rng(rng) { }

	errut::bool_t check(const Genome &g) override;
	errut::bool_t crossover(Genome &mutantDest, const Genome &origVector) override;
private:
	double m_CR;
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
errut::bool_t VectorDifferentialEvolutionMutation<T>::check(const Genome &g)
{
	const VectorGenome<T> *pGenome = dynamic_cast<const VectorGenome<T>*>(&g);
	if (!pGenome)
		return "Genome is not a VectorGenome of the correct type";
	return true;
}

template<class T>
inline std::shared_ptr<Genome> VectorDifferentialEvolutionMutation<T>::mutate(const Genome &r1, const Genome &r2, const Genome &r3)
{
	std::shared_ptr<Genome> result = r1.createCopy();
	VectorGenome<T> &g1 = static_cast<VectorGenome<T>&>(*result);
	const VectorGenome<T> &g2 = static_cast<const VectorGenome<T>&>(r2);
	const VectorGenome<T> &g3 = static_cast<const VectorGenome<T>&>(r3);
	
	std::vector<T> &v1 = g1.getValues();
	const std::vector<T> &v2 = g2.getValues();
	const std::vector<T> &v3 = g3.getValues();

	assert(v1.size() == v2.size());
	assert(v1.size() == v3.size());

	for (size_t i = 0 ; i < v1.size() ; i++)
		v1[i] += m_factor*(v2[i] - v3[i]);

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
inline errut::bool_t VectorDifferentialEvolutionCrossover<T>::crossover(Genome &mutantDest, const Genome &origVector)
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
		if (i != rndIdx && m_rng->getRandomDouble() > m_CR)
			v0[i] = v1[i];
	}
#else
	size_t j = (size_t)m_rng->getRandomUint32()%(v0.size());
	for (size_t k = 1 ; k <= v0.size() ; k++)
	{
		if (k != v0.size() && m_rng->getRandomDouble() > m_CR)
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

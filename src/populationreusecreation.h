#pragma once

#include "individual.h"
#include "population.h"
#include <list>

namespace eatk
{

// Idea is to let one EA continue with the population from another one
class PopulationReuseCreation : public eatk::IndividualCreation
{
public:
	PopulationReuseCreation(const std::shared_ptr<eatk::Population> &pop)
		: PopulationReuseCreation(std::vector<std::shared_ptr<eatk::Population>> { std::shared_ptr<eatk::Population>(pop) })
	{
	}

	PopulationReuseCreation(const std::vector<std::shared_ptr<eatk::Population>> &pops)
	{
		if (pops.size() == 0)
			return;

		for (auto &pop0 : pops)
		{
			const eatk::Population &pop = *pop0;
			if (pop.size() == 0)
				continue;

			m_referenceIndividual = pop.individual(0)->createCopy();
			break;
		}

		if (!m_referenceIndividual.get()) // Empty populations
			return; 

		m_referenceGenome = m_referenceIndividual->genomePtr()->createCopy();
		m_referenceFitness = m_referenceIndividual->fitnessPtr()->createCopy();

		for (auto &pop0 : pops)
		{
			const eatk::Population &pop = *pop0;
			for (auto &i : pop.individuals())
				m_genomePool.push_back(i->genomePtr()->createCopy());
		}

		//std::cerr << "DEBUG ReuseCreation: added " << m_genomePool.size() << " genomes to reuse" << std::endl;
	}

	std::shared_ptr<eatk::Genome> createInitializedGenome() override
	{
		if (m_genomePool.size() == 0)
			return nullptr;
		auto genome = m_genomePool.front();
		m_genomePool.pop_front();
		//std::cerr << "DEBUG ReuseCreation::createInitializedGenome called, " << m_genomePool.size() << " genomes remaining" << std::endl;
		return genome;
	}

	std::shared_ptr<eatk::Genome> createUnInitializedGenome() override
	{ 
		if (!m_referenceGenome.get())
			return nullptr;
		return m_referenceGenome->createCopy(false);
	}

	std::shared_ptr<eatk::Fitness> createEmptyFitness() override
	{
		if (!m_referenceFitness.get())
			return nullptr;
		return m_referenceFitness->createCopy(false);
	}

	std::shared_ptr<eatk::Individual> createReferenceIndividual() override { return m_referenceIndividual; }
private:
	std::list<std::shared_ptr<eatk::Genome>> m_genomePool;
	std::shared_ptr<eatk::Individual> m_referenceIndividual;
	std::shared_ptr<eatk::Genome> m_referenceGenome;
	std::shared_ptr<eatk::Fitness> m_referenceFitness;
};

} // end namespace

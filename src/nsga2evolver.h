#pragma once

#include "eatkconfig.h"
#include "populationevolver.h"
#include "crossovermutation.h"
#include "randomnumbergenerator.h"
#include "selection.h"
#include "singlepopulationcrossover.h"
#include "nondominatedsetcreator.h"
#include <memory>
#include <array>
#include <cmath>

namespace eatk
{

class NSGA2Evolver : public PopulationEvolver
{
public:
	NSGA2Evolver(
		const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<GenomeCrossover> &genomeCrossover,
		const std::shared_ptr<GenomeMutation> &genomeMutation,
		const std::shared_ptr<FitnessComparison> &fitComp, size_t numObjectives
		);
	~NSGA2Evolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_best; }
protected:
	virtual std::shared_ptr<NonDominatedSetCreator> allocatedNDSetCreator(const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives);
private:
	void buildWrapperPopulation(const Population &population);
	void calculateCrowdingDistances(const std::vector<std::shared_ptr<Individual>> &ndset) const;

	std::shared_ptr<FitnessComparison> m_fitComp;
	std::unique_ptr<SinglePopulationCrossover> m_crossover;
	std::shared_ptr<Population> m_tmpPop;
	const size_t m_numObjectives;

	class NSGA2IndividualWrapper : public Individual
	{
	public:
		NSGA2IndividualWrapper(size_t numObjectives,
		       const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			   size_t introducedInGeneration, size_t originalPosition)
			: Individual(genome, fitness, introducedInGeneration),
			  m_originalPosition(originalPosition)
		{
			m_fitnessDistances.resize(numObjectives);
		}

		std::shared_ptr<Individual> createNew(const std::shared_ptr<Genome> &genome, const std::shared_ptr<Fitness> &fitness,
			   size_t introducedInGeneration = std::numeric_limits<size_t>::max()) const override
		{
			return std::make_shared<NSGA2IndividualWrapper>(m_fitnessDistances.size(), genome, fitness, introducedInGeneration, std::numeric_limits<size_t>::max());
		}

		std::string toString() const override
		{
			std::string s;

			s += "pos: " + std::to_string(m_originalPosition) + " |";
			for (auto d : m_fitnessDistances)
				s += " " + std::to_string(d);
			s += "| " + Individual::toString();
			return s;
		}

		const size_t m_originalPosition;
		std::vector<double> m_fitnessDistances;
	};

	class NSGA2FitnessWrapper : public Fitness
	{
	public:
		NSGA2FitnessWrapper(const std::shared_ptr<Fitness> &origFitness) : m_origFitness(origFitness)
		{
			m_totalFitnessDistance = std::numeric_limits<double>::quiet_NaN();
			m_ndSetIndex = std::numeric_limits<size_t>::max();
		}

		std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override
		{
			std::shared_ptr<Fitness> realFitness = m_origFitness->createCopy(copyContents);
			std::shared_ptr<NSGA2FitnessWrapper> newFit = std::make_shared<NSGA2FitnessWrapper>(realFitness);
			if (copyContents)
			{
				newFit->m_totalFitnessDistance = m_totalFitnessDistance;
				newFit->m_ndSetIndex = m_ndSetIndex;
			}
			return newFit;
		}

		std::string toString() const override
		{
			return "{ ndset: " + std::to_string(m_ndSetIndex) + ", dist: " + std::to_string(m_totalFitnessDistance) + ", orig: " + m_origFitness->toString() + " }";
		}

		bool hasRealValues() const override { return m_origFitness->hasRealValues(); }
		double getRealValue(size_t objectiveNumber) const override { return m_origFitness->getRealValue(objectiveNumber); }

		std::shared_ptr<Fitness> m_origFitness;
		double m_totalFitnessDistance;
		size_t m_ndSetIndex;
	};

	class NSGA2FitnessWrapperOriginalComparison : public FitnessComparison
	{
	public:
		NSGA2FitnessWrapperOriginalComparison(const std::shared_ptr<FitnessComparison> &origCmp) : m_origCmp(origCmp) { }

		errut::bool_t check(const Fitness &f) const
		{
			if (!dynamic_cast<const NSGA2FitnessWrapper*>(&f))
				return "Expecting NSGA2FitnessWrapper object";
			const NSGA2FitnessWrapper &fw = static_cast<const NSGA2FitnessWrapper&>(f);
			return m_origCmp->check(*(fw.m_origFitness));
		}

		bool isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const
		{
			assert(dynamic_cast<const NSGA2FitnessWrapper*>(&first0) && dynamic_cast<const NSGA2FitnessWrapper*>(&second0));
			const NSGA2FitnessWrapper &first = static_cast<const NSGA2FitnessWrapper&>(first0);
			const NSGA2FitnessWrapper &second = static_cast<const NSGA2FitnessWrapper&>(second0);

			return m_origCmp->isFitterThan(*(first.m_origFitness), *(second.m_origFitness), objectiveNumber);
		}
	private:
		std::shared_ptr<FitnessComparison> m_origCmp;
	};

	class NSGA2FitWrapperNDSetCrowdingComparison : public FitnessComparison
	{
	public:
		errut::bool_t check(const Fitness &f) const
		{
			if (!dynamic_cast<const NSGA2FitnessWrapper*>(&f))
				return "Expecting NSGA2FitnessWrapper object";
			return true;
		}

		bool isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const
		{
			assert(dynamic_cast<const NSGA2FitnessWrapper*>(&first0) && dynamic_cast<const NSGA2FitnessWrapper*>(&second0));
			const NSGA2FitnessWrapper &first = static_cast<const NSGA2FitnessWrapper&>(first0);
			const NSGA2FitnessWrapper &second = static_cast<const NSGA2FitnessWrapper&>(second0);

			assert(first.m_ndSetIndex != std::numeric_limits<size_t>::max());
			assert(second.m_ndSetIndex != std::numeric_limits<size_t>::max());
			if (first.m_ndSetIndex < second.m_ndSetIndex)
				return true;
			if (second.m_ndSetIndex < first.m_ndSetIndex)
				return false;

			// Same ND set, check crowding distance
			assert(!std::isnan(first.m_totalFitnessDistance));
			assert(!std::isnan(second.m_totalFitnessDistance));

			// Use the less crowded one
			return first.m_totalFitnessDistance > second.m_totalFitnessDistance;
		}
	};

	std::vector<std::shared_ptr<Individual>> m_popWrapper;
	std::vector<std::shared_ptr<Individual>> m_best;
	std::shared_ptr<NSGA2FitnessWrapperOriginalComparison> m_fitOrigComp;
	std::shared_ptr<NonDominatedSetCreator> m_ndSetCreator;
};

}

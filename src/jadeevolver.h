#pragma once

#include "eatkconfig.h"
#include "differentialevolutionevolver.h"

namespace eatk
{

class JADEEvolver : public PopulationEvolver
{
public:
	JADEEvolver(const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber = 0, // negative means multi-objective
		double p = 0.05, double c = 0.1,
		bool useArchive = true,
		double initMuF = 0.5,
		double initMuCR = 0.5,
		size_t numObjectives = 1,
		const std::shared_ptr<NonDominatedSetCreator> &ndCreator = nullptr
		);
	~JADEEvolver();

	errut::bool_t check(const std::shared_ptr<Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize) override;
	const std::vector<std::shared_ptr<Individual>> &getBestIndividuals() const override { return m_bestIndividual; }
protected:
	virtual void onMutationCrossoverSettings(double muF, double muCR) const { }
private:
	void trimArchive(size_t targetPopulationSize);

	std::shared_ptr<RandomNumberGenerator> m_rng;
	std::shared_ptr<DifferentialEvolutionMutation> m_mut;
	std::shared_ptr<DifferentialEvolutionCrossover> m_cross;
	std::shared_ptr<FitnessComparison> m_fitComp;
	const int m_objectiveNumber;
	const size_t m_numObjectives;
	const double m_p, m_c;
	const double m_initMuF, m_initMuCR;
	const bool m_useArchive;
	double m_muF, m_muCR;
	std::vector<std::shared_ptr<Genome>> m_archive;

	std::vector<const Genome*> m_mutationGenomes;
	std::vector<double> m_mutationFactors;

	std::vector<double> m_CRi, m_Fi, m_SCR, m_SF;

	std::vector<std::shared_ptr<Individual>> m_bestIndividual;
	std::shared_ptr<NonDominatedSetCreator> m_ndCreator;
};

}

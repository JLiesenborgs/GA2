#pragma once

#include "population.h"

// TODO: MPI_Comm seems to be a pointer in OpenMPI

class MPIPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MPIPopulationFitnessCalculation();
	~MPIPopulationFitnessCalculation();

	// Layout will be exchanged between master and helpers, should be called
	// at same time at master and helpers
    // Need reference genome/fitness to be able to serialize it
	errut::bool_t init(const Genome &referenceGenome,
			    const Fitness &referenceFitness,
				std::shared_ptr<PopulationFitnessCalculation> &popCalc,
                MPI_Comm communicator,
                int root);

	errut::bool_t check(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t calculatePopulationFitness(const std::vector<std::shared_ptr<Population>> &populations) override;
	errut::bool_t calculatePopulationFitness_MPIHelper(); // should be called at same time as calculatePopulationFitness
private:
	std::shared_ptr<Genome> m_referenceGenome;
	std::shared_ptr<Fitness> m_referenceFitness;
	std::shared_ptr<Population> m_localPop;
	std::shared_ptr<PopulationFitnessCalculation> m_localPopulationFitnessCalculation;
    MPI_Comm m_comm;
    int m_root, m_mpiSize;

	std::vector<std::vector<std::pair<Genome *, Fitness *>>> m_helperGenomes;
};

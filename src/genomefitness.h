#pragma once

#include "mogal2config.h"
#include <mpi.h>
#include <errut/booltype.h>
#include <memory>

namespace mogal2
{

class Genome
{
public:
	Genome() { }
	virtual ~Genome() { }
	virtual std::shared_ptr<Genome> createCopy(bool copyContents = true) const { return nullptr; }
	virtual std::string toString() const { return "?"; };

	virtual errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) { return "Not implemented"; }
	virtual errut::bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const { return "Not implemented"; }
	virtual errut::bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) { return "Not implemented"; }
};

class Fitness
{
public:
	Fitness() : m_calculated(false) { }
	virtual ~Fitness() { }
	virtual std::shared_ptr<Fitness> createCopy(bool copyContents = true) const { return nullptr; }
	virtual std::string toString() const { return "?"; };
	bool isCalculated() const { return m_calculated; }
	void setCalculated(bool v = true) { m_calculated = v; }

	virtual errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) { return "Not implemented"; }
	virtual errut::bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const { return "Not implemented"; }
	virtual errut::bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) { return "Not implemented"; }
protected:
	bool m_calculated;
};

class FitnessComparison
{
public:
	FitnessComparison() { }
	virtual ~FitnessComparison() { }

	virtual errut::bool_t check(const Fitness &f) const { return "Not implemented in base class"; }
	// Using just a bool for speed - will be used in sorting functions
	// TODO: perhaps some kind of FitnessComparison operator is better?
	virtual bool isFitterThan(const Fitness &first, const Fitness &second, size_t objectiveNumber) const { return false; } // Implement this
};

class GenomeFitnessCalculation
{
public:
	GenomeFitnessCalculation() { }
	virtual ~GenomeFitnessCalculation() { }

	// These are to allow a more async version, but by default the sync version is called
	virtual errut::bool_t startNewCalculation(const Genome &genome) { return true; }
	
	// Return type says if something went wrong; fitness.isCalculated marks when done
	// May need multiple calls, to make async behaviour possible
	virtual errut::bool_t pollCalculate(const Genome &genome, Fitness &fitness)
	{
		auto r = calculate(genome, fitness);
		fitness.setCalculated();
		return r;
	}

	// Convenience function for sync operation
	virtual errut::bool_t calculate(const Genome &genome, Fitness &fitness) { return "Not implemented"; }
};

class GenomeFitnessCreation
{
public:
	GenomeFitnessCreation() { }
	virtual ~GenomeFitnessCreation() { }

    // TODO: how to signal error?
    virtual std::shared_ptr<Genome> createInitializedGenome() = 0;
    virtual std::shared_ptr<Fitness> createEmptyFitness() = 0;
};

}
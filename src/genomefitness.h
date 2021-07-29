#pragma once

#include "eatkconfig.h"

#ifdef EATKCONFIG_MPISUPPORT
#include <mpi.h>
#endif // EATKCONFIG_MPISUPPORT

#include <errut/booltype.h>
#include <memory>
#include <vector>

namespace eatk
{

template<class T>
class GenomeFitnessBase
{
protected:
	GenomeFitnessBase() { }
public:
	virtual ~GenomeFitnessBase() { }

	virtual std::string toString() const { return "?"; };

	virtual std::shared_ptr<T> createCopy(bool copyContents = true) const { return nullptr; }

#ifdef EATKCONFIG_MPISUPPORT
	virtual errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) { return "Not implemented"; }
	virtual errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator,
								   std::vector<MPI_Request> &requests) const { return "Not implemented"; }
	virtual errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator,
								   std::vector<MPI_Request> &requests) { return "Not implemented"; }
#endif // EATKCONFIG_MPISUPPORT
};

class Genome : public GenomeFitnessBase<Genome>
{
public:
	Genome() { }
	~Genome() { }
};

class Fitness : public GenomeFitnessBase<Fitness>
{
public:
	Fitness() : m_calculated(false) { }
	~Fitness() { }

	bool isCalculated() const { return m_calculated; }
	void setCalculated(bool v = true) { m_calculated = v; }
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

	// TODO: do we need a check function here?

	// genomesForPopulationCalculator is e.g. the number that needs to be calculated across threads
	// The MPI implementation itself uses a different local calculator, so doesn't call this itself
	// This means that it's mainly a multi-thread thing
	virtual errut::bool_t onNewCalculationStart(size_t genomesForThisCalculator, size_t genomesForPopulationCalculator)  { return true; }
	virtual errut::bool_t onCalculationStarted() { return true; }
	virtual errut::bool_t onCalculationEnded() { return true; }

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

}
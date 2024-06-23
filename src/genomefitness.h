#pragma once

#include "eatkconfig.h"

#ifdef EATKCONFIG_MPISUPPORT
#include <mpi.h>
#endif // EATKCONFIG_MPISUPPORT

#include <errut/booltype.h>
#include <memory>
#include <vector>
#include <functional>
#include <limits>

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

	// TODO: is float better?
	// This is needed in case the working of an algorithm requires the
	// fitness value to be numerical (NSGA2 needs this for example)
	virtual bool hasRealValues() const { return false; }
	virtual double getRealValue(size_t objectiveNumber) const { return std::numeric_limits<double>::quiet_NaN(); }

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

	static std::function<bool(const Fitness &f1, const Fitness &f2)> getDominatesFunction(const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber, size_t numObjectives);
};

inline std::function<bool(const Fitness &f1, const Fitness &f2)> FitnessComparison::getDominatesFunction(const std::shared_ptr<FitnessComparison> &fitComp, int objectiveNumber, size_t numObjectives)
{
	std::function<bool(const Fitness &f1, const Fitness &f2)> isFitterThan;

	if (numObjectives <= 1) // single objective
	{
		isFitterThan = [fitComp, objectiveNumber](const Fitness &f1, const Fitness &f2)
		{
			return fitComp->isFitterThan(f1, f2, objectiveNumber);
		};
	}
	else
	{
		// multi-objective, use dominance
		isFitterThan = [fitComp, numObjectives](const Fitness &f1, const Fitness &f2)
		{
			size_t betterOrEqualCount = 0;
			size_t betterCount = 0;
			for (size_t i = 0 ; i < numObjectives ; i++)
			{
				if (fitComp->isFitterThan(f1, f2, i))
				{
					betterCount++;
					betterOrEqualCount++;
				}
				else // f1 not strictly better than f2 for i
				{
					if (!fitComp->isFitterThan(f2, f1, i)) // then they must have equal fitness
					{
						betterOrEqualCount++;
					}
					else
					{
						// We can never get betterOrEqualCount == m_numObjectives
						return false;
					}
				}
			}
			// if we got here, then betterOrEqualCount == m_numObjectives
			if (betterCount > 0)
				return true;
			return false;
		};
	}
	return isFitterThan;
}

}

#pragma once

#include "eatkconfig.h"
#include "genomefitness.h"
#include <memory>
#include <string>
#include <sstream>

namespace eatk
{

template<class T>
class ValueFitness : public Fitness
{
public:
	ValueFitness() { }
	ValueFitness(T value) : m_value(value) { }
	~ValueFitness() { }

	T getValue() const { return m_value; }
	void setValue(T x) { m_value = x; }

	bool hasRealValues() const override { return true; }
	double getRealValue(size_t objectiveNumber) const override { return (double)getValue(); }

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		auto f = std::make_shared<ValueFitness<T>>(m_value);
		if (copyContents && Fitness::isCalculated())
			f->setCalculated();
		return f;
	}

	std::string toString() const override
	{
		if (!Fitness::isCalculated())
			return Fitness::toString();
#if 1
		return std::to_string(m_value);
#else
		std::stringstream ss;
		ss.precision(15);
		ss << m_value;
		return ss.str();
#endif
	}

#ifdef EATKCONFIG_MPISUPPORT
	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override
	{
		// Nothing to do here - I think we can safely omit this
		return true;
	}

	errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator, std::vector<MPI_Request> &requests) const override
	{
		requests.resize(1);
		MPI_Isend(&m_value, 1, m_mpiType, dest, tag, communicator, requests.data());
		// ::MPI_Send(&m_value, 1, m_mpiType, dest, tag, communicator);
		return true;
	}

	errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator, std::vector<MPI_Request> &requests) override
	{
		requests.resize(1);
		MPI_Irecv(&m_value, 1, m_mpiType, src, tag, communicator, requests.data());
		// ::MPI_Recv(&m_value, 1, m_mpiType, src, tag, communicator, MPI_STATUS_IGNORE);
		return true;
	}
#endif // EATKCONFIG_MPISUPPORT
private:
	T m_value;
#ifdef EATKCONFIG_MPISUPPORT
	static MPI_Datatype m_mpiType;
#endif // EATKCONFIG_MPISUPPORT
};

template<class T, bool minimum = true>
class ValueFitnessComparison : public FitnessComparison
{
public:
	ValueFitnessComparison() { }
	~ValueFitnessComparison() { }

	errut::bool_t check(const Fitness &f) const override
	{
		const ValueFitness<T> *pVf = dynamic_cast<const ValueFitness<T>*>(&f);
		if (pVf)
		 	return true;

		return "Different fitness type than expected: " + std::string(typeid(f).name());
	}

	bool isFitterThan(const Fitness &first, const Fitness &second, size_t objectiveNumber) const override
	{
		const ValueFitness<T> &f1 = static_cast<const ValueFitness<T> &>(first);
		const ValueFitness<T> &f2 = static_cast<const ValueFitness<T> &>(second);

		if (minimum)
			return (f1.getValue() < f2.getValue());
		return (f1.getValue() > f2.getValue());
	}
};

}

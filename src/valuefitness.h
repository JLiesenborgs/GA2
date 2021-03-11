#pragma once

#include "mogal2config.h"
#include "genomefitness.h"
#include <memory>
#include <string>

namespace mogal2
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
        return std::to_string(m_value);
    }

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override
    {
        // Nothing to do here - I think we can safely omit this
        return true;
    }

	errut::bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const override
    {
        MPI_Isend(&m_value, 1, m_mpiType, dest, tag, communicator, pRequest);
        return true;
    }

	errut::bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) override
    {
        MPI_Irecv(&m_value, 1, m_mpiType, src, tag, communicator, pRequest);
        return true;
    }
private:
    T m_value;
    static MPI_Datatype m_mpiType;
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

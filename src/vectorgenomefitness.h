#pragma once

#include "eatkconfig.h"
#include "valuevector.h"
#include "genomefitness.h"
#include <cassert>

namespace eatk
{

template<class T>
class VectorGenome : public ValueVector<Genome, T>
{
public:
	VectorGenome(size_t n = 0) : ValueVector<Genome,T>(n) { }
	VectorGenome(size_t n, T initValue) : ValueVector<Genome,T>(n, initValue) { }
	~VectorGenome() { }

	std::shared_ptr<Genome> createCopy(bool copyContents = true) const override
	{
		return ValueVector<Genome, T>::template createCopy<VectorGenome<T>>(copyContents);
	}

// For DE-like crossover/mutation in GA
	static void setValue(Genome &f, size_t pos, T value)
	{
		static_cast<VectorGenome<T>&>(f).m_values[pos] = value;
	}

	static T getValue(const Genome &f, size_t pos)
	{
		return static_cast<const VectorGenome<T>&>(f).m_values[pos];
	}

	static size_t getSize(const Genome &f)
	{
		return static_cast<const VectorGenome<T>&>(f).m_values.size();
	}
};

template<class T>
class VectorFitness : public ValueVector<Fitness, T>
{
public:
	VectorFitness(size_t n = 0) : ValueVector<Fitness, T>(n, (T)0) { }
	~VectorFitness() { }

	std::string toString() const override
	{
		if (!Fitness::isCalculated())
			return "?";
		return ValueVector<Fitness, T>::toString();
	}

	bool hasRealValues() const override { return true; }
	double getRealValue(size_t objectiveNumber) const override
	{
		assert(objectiveNumber < (ValueVector<Fitness, T>::getValues().size()));
		return (double)ValueVector<Fitness, T>::getValues()[objectiveNumber];
	}

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		auto g = ValueVector<Fitness, T>::template createCopy<VectorFitness<T>>(copyContents);
		if (copyContents && Fitness::isCalculated())
			g->setCalculated();
		return g;
	}
};

template<class T, bool minimum = true>
class VectorFitnessComparison : public FitnessComparison
{
public:
	VectorFitnessComparison() { }
	~VectorFitnessComparison() { }

	errut::bool_t check(const Fitness &f) const override
	{
		const VectorFitness<T> *pVf = dynamic_cast<const VectorFitness<T>*>(&f);
		if (pVf)
		 	return true;

		return "Different fitness type than expected: " + std::string(typeid(f).name());
	}

	bool isFitterThan(const Fitness &first, const Fitness &second, size_t objectiveNumber) const override
	{
		const VectorFitness<T> &f1 = static_cast<const VectorFitness<T> &>(first);
		const VectorFitness<T> &f2 = static_cast<const VectorFitness<T> &>(second);

		const std::vector<T> &v1 = f1.getValues();
		const std::vector<T> &v2 = f2.getValues();
		assert(objectiveNumber >= 0 && objectiveNumber < v1.size());
		assert(v1.size() == v2.size());

		if constexpr (minimum)
			return (v1[objectiveNumber] < v2[objectiveNumber]);
		else
			return (v1[objectiveNumber] > v2[objectiveNumber]);
	}
};

typedef VectorGenome<float> FloatVectorGenome;
typedef VectorFitness<float> FloatVectorFitness;

typedef VectorGenome<double> DoubleVectorGenome;
typedef VectorFitness<double> DoubleVectorFitness;

typedef VectorGenome<int> IntVectorGenome;
typedef VectorFitness<int> IntVectorFitness;

}

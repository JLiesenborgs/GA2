#pragma once

#include "valuevector.h"
#include "genomefitness.h"
#include <cassert>

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

		if (minimum)
			return (v1[objectiveNumber] < v2[objectiveNumber]);
		return (v1[objectiveNumber] > v2[objectiveNumber]);
	}
};

typedef VectorGenome<float> FloatVectorGenome;
typedef VectorFitness<float> FloatVectorFitness;

typedef VectorGenome<double> DoubleVectorGenome;
typedef VectorFitness<double> DoubleVectorFitness;

typedef VectorGenome<int> IntVectorGenome;
typedef VectorFitness<int> IntVectorFitness;


#pragma once

#include "valuevector.h"
#include "genomefitness.h"
#include <cassert>

template<class T>
class VectorGenome : public ValueVector<Genome, float>
{
public:
	VectorGenome(size_t n = 0) : ValueVector(n) { }
	VectorGenome(size_t n, float initValue) : ValueVector(n, initValue) { }
	~VectorGenome() { }

	std::shared_ptr<Genome> createCopy(bool copyContents = true) const override
	{
		return ValueVector::createCopy<VectorGenome<T>>(copyContents);
	}
};

template<class T>
class VectorFitness : public ValueVector<Fitness, float>
{
public:
	VectorFitness(size_t n = 0) : ValueVector(n, 0) { }
	~VectorFitness() { }

	bool isFitterThan(const Fitness &otherFitness, int objectiveNumber) const override
	{
		assert(objectiveNumber >= 0 && objectiveNumber < (int)m_values.size());
		const VectorFitness<T> &otherVectorFitness = static_cast<const VectorFitness<T> &>(otherFitness);
		
		if (m_values[objectiveNumber] < otherVectorFitness.m_values[objectiveNumber])
			return true;
		return false;
	}

	std::string toString() const override
	{
		if (!isCalculated())
			return "?";
		return ValueVector::toString();
	}

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		auto g = ValueVector::createCopy<VectorFitness<T>>(copyContents);
		if (copyContents && isCalculated())
			g->setCalculated();
		return g;
	}
};

typedef VectorGenome<float> FloatVectorGenome;
typedef VectorFitness<float> FloatVectorFitness;

typedef VectorGenome<double> DoubleVectorGenome;
typedef VectorFitness<double> DoubleVectorFitness;

typedef VectorGenome<int> IntVectorGenome;
typedef VectorFitness<int> IntVectorFitness;


#pragma once

#include "valuevector.h"
#include "genomefitness.h"

class FloatVectorGenome : public ValueVector<Genome, float>
{
public:
	FloatVectorGenome(size_t n = 0);
	FloatVectorGenome(size_t n, float initValue);
	~FloatVectorGenome();

	std::shared_ptr<Genome> createCopy(bool copyContents = true) const override;
};

class FloatVectorFitness : public ValueVector<Fitness, float>
{
public:
	FloatVectorFitness(size_t n = 0);
	~FloatVectorFitness();

	std::string toString() const override;
	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override;
};

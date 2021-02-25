#include "floatvectorgenomefitness.h"

using namespace std;

template<> MPI_Datatype ValueVector<Genome, float>::m_mpiType = MPI_FLOAT;
template<> MPI_Datatype ValueVector<Fitness, float>::m_mpiType = MPI_FLOAT;

FloatVectorGenome::FloatVectorGenome(size_t n) : ValueVector(n)
{ 
}

FloatVectorGenome::FloatVectorGenome(size_t n, float initValue) : ValueVector(n, initValue)
{
}

FloatVectorGenome::~FloatVectorGenome()
{
}

shared_ptr<Genome> FloatVectorGenome::createCopy(bool copyContents) const
{
    return ValueVector::createCopy<FloatVectorGenome>(copyContents);
}	

FloatVectorFitness::FloatVectorFitness(size_t n) : ValueVector(n, 0)
{
}

FloatVectorFitness::~FloatVectorFitness()
{
}

string FloatVectorFitness::toString() const
{
    if (!isCalculated())
        return "?";
    return ValueVector::toString();
}

shared_ptr<Fitness> FloatVectorFitness::createCopy(bool copyContents) const
{
    auto g = ValueVector::createCopy<FloatVectorFitness>(copyContents);
    if (copyContents && isCalculated())
        g->setCalculated();
    return g;
}

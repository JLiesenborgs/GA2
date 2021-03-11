#include "vectorgenomefitness.h"

using namespace std;

namespace mogal2
{

template<> MPI_Datatype ValueVector<Genome, float>::m_mpiType = MPI_FLOAT;
template<> MPI_Datatype ValueVector<Fitness, float>::m_mpiType = MPI_FLOAT;

template<> MPI_Datatype ValueVector<Genome, double>::m_mpiType = MPI_DOUBLE;
template<> MPI_Datatype ValueVector<Fitness, double>::m_mpiType = MPI_DOUBLE;

template<> MPI_Datatype ValueVector<Genome, int>::m_mpiType = MPI_INT;
template<> MPI_Datatype ValueVector<Fitness, int>::m_mpiType = MPI_INT;

}

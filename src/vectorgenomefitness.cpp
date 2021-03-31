#include "vectorgenomefitness.h"

using namespace std;

namespace eatk
{

template<> MPI_Datatype ValueVector<Genome, float>::m_mpiType = MPI_FLOAT;
template<> MPI_Datatype ValueVector<Fitness, float>::m_mpiType = MPI_FLOAT;

template<> MPI_Datatype ValueVector<Genome, double>::m_mpiType = MPI_DOUBLE;
template<> MPI_Datatype ValueVector<Fitness, double>::m_mpiType = MPI_DOUBLE;

template<> MPI_Datatype ValueVector<Genome, int>::m_mpiType = MPI_INT;
template<> MPI_Datatype ValueVector<Fitness, int>::m_mpiType = MPI_INT;

template<> MPI_Datatype ValueVector<Genome, unsigned int>::m_mpiType = MPI_UNSIGNED;
template<> MPI_Datatype ValueVector<Fitness, unsigned int>::m_mpiType = MPI_UNSIGNED;

template<> MPI_Datatype ValueVector<Genome, char>::m_mpiType = MPI_CHAR;
template<> MPI_Datatype ValueVector<Fitness, char>::m_mpiType = MPI_CHAR;

template<> MPI_Datatype ValueVector<Genome, unsigned char>::m_mpiType = MPI_UNSIGNED_CHAR;
template<> MPI_Datatype ValueVector<Fitness, unsigned char>::m_mpiType = MPI_UNSIGNED_CHAR;

}

#include "valuefitness.h"

#ifdef EATKCONFIG_MPISUPPORT

namespace eatk
{

template<> MPI_Datatype ValueFitness<float>::m_mpiType = MPI_FLOAT;
template<> MPI_Datatype ValueFitness<double>::m_mpiType = MPI_DOUBLE;
template<> MPI_Datatype ValueFitness<int>::m_mpiType = MPI_INT;
template<> MPI_Datatype ValueFitness<unsigned int>::m_mpiType = MPI_UNSIGNED;
template<> MPI_Datatype ValueFitness<char>::m_mpiType = MPI_CHAR;
template<> MPI_Datatype ValueFitness<unsigned char>::m_mpiType = MPI_UNSIGNED_CHAR;

}

#endif // EATKCONFIG_MPISUPPORT
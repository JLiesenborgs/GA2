#include "valuefitness.h"

template<> MPI_Datatype ValueFitness<float>::m_mpiType = MPI_FLOAT;
template<> MPI_Datatype ValueFitness<double>::m_mpiType = MPI_DOUBLE;
template<> MPI_Datatype ValueFitness<int>::m_mpiType = MPI_INT;

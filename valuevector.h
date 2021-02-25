#pragma once

#include <mpi.h>
#include <errut/booltype.h>
#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <cassert>
#include <sstream>

template<class Base, class Type>
class ValueVector : public Base
{
public:
	ValueVector(size_t n = 0) : m_values(n) { }
	ValueVector(size_t n, Type initValue) : m_values(n, initValue) { }
	~ValueVector() { }

	std::vector<Type> &getValues() { return m_values; }
	const std::vector<Type> &getValues() const { return m_values; }

	std::string toString() const override
	{
		std::stringstream ss;

		ss << "[";
		for (auto x : m_values)
			ss << " " << x;
		ss << " ]";
		return ss.str();
	}

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override
	{
		int num = m_values.size();
		MPI_Bcast(&num, 1, MPI_INT, root, communicator);
		m_values.resize(num);
		return true;
	}

	errut::bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const override
	{
		// TODO: does this work when m_values.size() == 0? Should it?

		// Master and helper should already know the genome layout, no need to send the
		// number of values first
		MPI_Isend(m_values.data(), m_values.size(), m_mpiType, dest, tag, communicator, pRequest);
		return true;
	}

	errut::bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) override
	{
		// cerr << "Receiving " << m_values.size() << " floats" << endl;
		MPI_Irecv(m_values.data(), m_values.size(), m_mpiType, src, tag, communicator, pRequest);
		return true;
	}

	template<class Derived>
	std::shared_ptr<Derived> createCopy(bool copyContents = true) const
	{
		auto g = std::make_shared<Derived>(m_values.size());
		if (copyContents && m_values.size() > 0)
		{
			assert(m_values.size() == g->m_values.size());
			std::memcpy(g->m_values.data(), m_values.data(), sizeof(Type)*m_values.size());
		}
		return g;
	}
protected:
	std::vector<Type> m_values;
	static MPI_Datatype m_mpiType;
};

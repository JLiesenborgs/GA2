#pragma once

#include "eatkconfig.h"
#include <errut/booltype.h>
#include <mpi.h>
#include <memory>
#include <vector>

namespace eatk
{

// To avoid circular references, no shared_ptr is used, but a weak_ptr

class MPIEventHandler
{
public:
    enum EventType { Done, Calculation, MaxEvent };

    MPIEventHandler() { }
    virtual ~MPIEventHandler() { }
    virtual errut::bool_t handleEvent(EventType t) { return "Not implemented in base class"; }
};

class MPIEventDistributor
{
public:
    MPIEventDistributor(int root = 0, MPI_Comm communicator = MPI_COMM_WORLD);
    ~MPIEventDistributor();

    errut::bool_t setHandler(MPIEventHandler::EventType evtType, const std::weak_ptr<MPIEventHandler> &calc);
    errut::bool_t signal(MPIEventHandler::EventType evtType);
    errut::bool_t eventLoop(); // Should be started on the non-root members of the communicator
private:
    MPI_Comm m_comm;
    int m_root;
    std::vector<std::weak_ptr<MPIEventHandler>> m_handlers;
};

}
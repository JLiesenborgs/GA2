#include "mpieventdistributor.h"
#include <string>

using namespace errut;
using namespace std;

namespace eatk
{

MPIEventDistributor::MPIEventDistributor(int root, MPI_Comm communicator)
    : m_comm(communicator), m_root(root)
{
    m_handlers.resize(MPIEventHandler::MaxEvent);
}

MPIEventDistributor::~MPIEventDistributor()
{
}

errut::bool_t MPIEventDistributor::setHandler(MPIEventHandler::EventType evtType,
                                              const std::weak_ptr<MPIEventHandler> &calc)
{
    int evt = (int)evtType;
    if (evt < 0 || evt >= (int)m_handlers.size())
        return "Invalid event type " + to_string(evt);
    m_handlers[evt] = calc;
    return true;
}

bool_t MPIEventDistributor::signal(MPIEventHandler::EventType evtType)
{
    int evt = (int)evtType;
    MPI_Bcast(&evt, 1, MPI_INT, m_root, m_comm);
    return true;
}

bool_t MPIEventDistributor::eventLoop()
{
    while (true)
    {
        int evt = -1; // Set this to something invalid, in case we're misidentifying ourselves as root
        MPI_Bcast(&evt, 1, MPI_INT, m_root, m_comm);

        if ((MPIEventHandler::EventType)evt == MPIEventHandler::Done)
            break;

        if (evt < 0 || evt >= (int)m_handlers.size())
            return "Received invalid event type " + to_string(evt);

        if (auto handler = m_handlers[evt].lock())
        {
            bool_t r = handler->handleEvent((MPIEventHandler::EventType)evt);       
            if (!r)
                return "Error handling event " + to_string(evt) + ": " + r.getErrorString();
        }
        else
            return "No handler set for event type " + to_string(evt) + "(or handler expired)";
    }
    return true;
}

}

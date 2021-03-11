#include "mpipopulationfitnesscalculation.h"
#include <iostream>
#include <cassert>
#include <typeinfo>

using namespace errut;
using namespace std;

namespace mogal2
{

MPIPopulationFitnessCalculation::MPIPopulationFitnessCalculation(const std::weak_ptr<MPIEventDistributor> &mpiDist)
	: m_evtDist(mpiDist)
{
}

MPIPopulationFitnessCalculation::~MPIPopulationFitnessCalculation()
{
}

bool_t MPIPopulationFitnessCalculation::init(const Genome &referenceGenome, const Fitness &referenceFitness,
											 shared_ptr<PopulationFitnessCalculation> popCalc,
											 MPI_Comm communicator, int root)
{
	int size = 0;
	MPI_Comm_size(communicator, &size);
	if (root < 0 || root >= size)
		return "Invalid root (" + to_string(root) + "), should be less than " + to_string(size);
	
	m_mpiSize = size;
	m_root = root;
	m_comm = communicator;

	m_referenceGenome = referenceGenome.createCopy(false);
	m_referenceFitness = referenceFitness.createCopy(false);
	m_localPopulationFitnessCalculation = popCalc;
	m_localPop = make_shared<Population>();

	bool_t r = m_referenceGenome->MPI_BroadcastLayout(m_root, m_comm);
	if (!r)
		return "Error broadcasting genome layout: " + r.getErrorString();
	if (!(r = m_referenceFitness->MPI_BroadcastLayout(m_root, m_comm)))
		return "Error broadcasting fitness layout: " + r.getErrorString();

	return true;
}

bool_t MPIPopulationFitnessCalculation::check(const std::vector<std::shared_ptr<Population>> &populations)
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";
	if (!m_localPopulationFitnessCalculation.get())
		return "Local fitness calculation not set";
	
	for (auto &pop : populations)
	{
		for (auto &i : pop->individuals())
		{
			string popGenomeType = typeid(i->genomeRef()).name();
			string refGenomeType = typeid(*(m_referenceGenome.get())).name();
			if (popGenomeType != refGenomeType)
				return "Genome in population is of different type than reference genome (" + popGenomeType + " != " + refGenomeType +  ")";

			string popFitnessType = typeid(i->fitnessRef()).name();
			string refFitnessType = typeid(*(m_referenceFitness.get())).name();
			if (popFitnessType != refFitnessType)
				return "Fitness in population is of different type than reference fitness (" + popFitnessType + " != " + refFitnessType +  ")";
		}
	}
	return true;
}

// TODO: send/check genome/fitness type, to make sure helper is using correct types
bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";
	if (!m_localPopulationFitnessCalculation.get())
		return "Local fitness calculation not set";

	if (auto dist = m_evtDist.lock())
		dist->signal(MPIEventHandler::Calculation);

	// TODO: check reference fitness type/layout against population?

	m_localPop->clear();

	m_helperGenomes.resize(m_mpiSize); // TODO: use a reusable array
	for (auto &h : m_helperGenomes)
		h.clear();

	int nextHelper = 0;
	int localCount = 0;
	int remoteCount = 0;

	// Split the work over the helpers
	for (auto pPop : populations)
	{		
		for (auto &i : pPop->individuals())
		{
			// TODO: check for first that type is same as reference genome? Or only for debugging?
			// cerr << typeid(*(i->m_genome.get())).name() << " " << typeid(*(m_referenceGenome.get())).name() << endl;

			if (!i->fitnessRef().isCalculated())
			{
				if (nextHelper == m_root) // This is for local calculation
				{
					m_localPop->append(i);
					localCount++;
				}
				else
				{
					m_helperGenomes[nextHelper].push_back({i->genomePtr(), i->fitnessPtr()});
					remoteCount++;
				}

				nextHelper = (nextHelper + 1)%m_mpiSize;
			}
		}
	}

	vector<MPI_Request> requests;
	auto getNextRequest = [&requests]()
	{
		requests.push_back(MPI_Request());
		return &(requests.back());
	};

	// First send everything to the helpers
	for (int helper = 0 ; helper < (int)m_helperGenomes.size() ; helper++)
	{
		if (helper == m_root) // Don't send the stuff we're doing locally
			continue;
		
		auto &thisHelperGenomes = m_helperGenomes[helper];
		int numGenomes = thisHelperGenomes.size();
		MPI_Send(&numGenomes, 1, MPI_INT, helper, 0, m_comm);

		for (int individual = 0 ; individual < numGenomes ; individual++)
		{
			Genome *pGenome = thisHelperGenomes[individual].first;
			MPI_Request *pReq = getNextRequest();
			auto r = pGenome->MPI_ISend(helper, individual, m_comm, pReq);

			if (!r)
				return "Error in genome's MPI_ISend: " + r.getErrorString();
		}
	}

	// Wait for the work to be received (we may not be doing MPI calls for a while,
	// and don't want the sending to stall)
	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	
	// Do our own work
	auto r = m_localPopulationFitnessCalculation->calculatePopulationFitness({ m_localPop });
	if (!r)
		return "Error in local fitness calculation: " + r.getErrorString();

	// TODO: check calculation flags in localpop?
	
	// Receive the calculations
	requests.clear();
	for (int helper = 0 ; helper < (int)m_helperGenomes.size() ; helper++)
	{
		auto &partHelperGenoms = m_helperGenomes[helper];
		for (int individual = 0 ; individual < (int)partHelperGenoms.size() ; individual++)
		{
			Fitness *pFitness = partHelperGenoms[individual].second;
			assert(pFitness);

			pFitness->MPI_IRecv(helper, individual, m_comm, getNextRequest());
			pFitness->setCalculated();
			// cerr << "Starting receive for " << helper << "," << individual << endl;
		}
	}

	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	//cerr << "Local: " << localCount << ", remote: " << remoteCount << endl;

	return true;
}

bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness_MPIHelper()
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";
	if (!m_localPopulationFitnessCalculation.get())
		return "Local fitness calculation not set";

	int numGenomes = 0;
	MPI_Recv(&numGenomes, 1, MPI_INT, m_root, 0, m_comm, MPI_STATUS_IGNORE);
	// cerr << "Calculating " << numGenomes << " in helper" << endl;

	// TODO: recycle previously created instances
	vector<pair<Genome *, Fitness *>> helperGenomes;

	m_localPop->clear();
	for (int i = 0 ; i < numGenomes ; i++) 
	{
		auto genome = m_referenceGenome->createCopy(false);
		auto fitness = m_referenceFitness->createCopy(false);
		fitness->setCalculated(false); // Make sure it's not set accidentally, so that we actually calculate it

		m_localPop->append(make_shared<Individual>(genome, fitness));
		helperGenomes.push_back({genome.get(), fitness.get()});
	}

	vector<MPI_Request> requests(numGenomes);
	for (int i = 0 ; i < numGenomes ; i++)
	{
		Genome *pGenome = helperGenomes[i].first;
		auto r = pGenome->MPI_IRecv(m_root, i, m_comm, &requests[i]); 
		if (!r)
			return "Error receiving genome in helper: " + r.getErrorString();
	}

	// Wait for everything to be received
	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	// cerr << "Received genomes in helper" << endl;
	
	// TODO: make sure that fitness calculation flags are cleared!!
	auto r = m_localPopulationFitnessCalculation->calculatePopulationFitness({ m_localPop });
	if (!r)
		return "Error in local fitness calculation: " + r.getErrorString();

	// Send fitness results back
	for (int individual = 0 ; individual < (int)helperGenomes.size() ; individual++)
	{
		Fitness *pFitness = helperGenomes[individual].second;
		if (!(r = pFitness->MPI_ISend(m_root, individual, m_comm, &requests[individual])))
			return "Error sending back fitness: " + r.getErrorString();
		
		// cerr << "Sending back fitness for " << individual << ": " << pFitness->toString() << endl;

		pFitness->setCalculated(false); // Already clear flag again
	}

	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	
	return true;
}

bool_t MPIPopulationFitnessCalculation::handleEvent(MPIEventHandler::EventType t)
{
	if (t != MPIEventHandler::Calculation)
		return "Can't handle event " + to_string(t);
	return calculatePopulationFitness_MPIHelper();
}

}

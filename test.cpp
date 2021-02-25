#include "population.h"
#include <vector>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <cstring>
#include <cassert>
#include <cstdio>

// TODO: namespace

// TODO: we're assuming that rank 0 is the master

// TODO: make communicator configurable

// TODO: for genome and fitness, a copyInto function? So that memory doesn't need to be
//       allocated/reallocated often

// TODO: MPI_Comm seems to be a pointer in OpenMPI

using namespace std;
using namespace errut;

template<class Base>
class FloatVector : public Base
{
public:
	FloatVector(size_t n = 0) : m_values(n) { }
	FloatVector(size_t n, float initValue) : m_values(n, initValue) { }
	~FloatVector() { }

	vector<float> &getValues() { return m_values; }
	const vector<float> &getValues() const { return m_values; }

	string toString() const override
	{
		stringstream ss;

		ss << "[";
		for (auto x : m_values)
			ss << " " << x;
		ss << " ]";
		return ss.str();
	}

	bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override
	{
		int num = m_values.size();
		MPI_Bcast(&num, 1, MPI_INT, root, communicator);
		m_values.resize(num);
		return true;
	}

	bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const override
	{
		// TODO: does this work when m_values.size() == 0? Should it?

		// Master and helper should already know the genome layout, no need to send the
		// number of values first
		MPI_Isend(m_values.data(), m_values.size(), MPI_FLOAT, dest, tag, communicator, pRequest);
		return true;
	}

	bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) override
	{
		// cerr << "Receiving " << m_values.size() << " floats" << endl;
		MPI_Irecv(m_values.data(), m_values.size(), MPI_FLOAT, src, tag, communicator, pRequest);
		return true;
	}

	template<class Derived>
	shared_ptr<Derived> createCopy(bool copyContents = true) const
	{
		auto g = make_shared<Derived>(m_values.size());
		if (copyContents && m_values.size() > 0)
		{
			assert(m_values.size() == g->m_values.size());
			memcpy(g->m_values.data(), m_values.data(), sizeof(float)*m_values.size());
		}
		return g;
	}
protected:
	vector<float> m_values;
};

class FloatVectorGenome : public FloatVector<Genome>
{
public:
	FloatVectorGenome(size_t n = 0) : FloatVector<Genome>(n) { }
	FloatVectorGenome(size_t n, float initValue) : FloatVector<Genome>(n, initValue) { }
	~FloatVectorGenome() { }

	shared_ptr<Genome> createCopy(bool copyContents = true) const override
	{
		return FloatVector<Genome>::createCopy<FloatVectorGenome>(copyContents);
	}	
};

class FloatVectorFitness : public FloatVector<Fitness>
{
public:
	FloatVectorFitness(size_t n = 0) : FloatVector<Fitness>(n, 0) { }
	~FloatVectorFitness() { }

	string toString() const override
	{
		if (!isCalculated())
			return "?";
		return FloatVector<Fitness>::toString();
	}

	shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		auto g = FloatVector<Fitness>::createCopy<FloatVectorFitness>(copyContents);
		if (copyContents && isCalculated())
			g->setCalculated();
		return g;
	}
};

class DummyGenomeFitnessCalculation : public GenomeFitnessCalculation
{
public:
	DummyGenomeFitnessCalculation() { }
	~DummyGenomeFitnessCalculation() { }
	bool_t pollCalculate(const Genome &genome, Fitness &fitness)
	{
		// TODO: move checking code to separate routine
		const FloatVectorGenome *pGenome = dynamic_cast<const FloatVectorGenome *>(&genome);
		if (!pGenome)
			return "Genome is not of expected type";
		FloatVectorFitness *pFitness = dynamic_cast<FloatVectorFitness *>(&fitness);
		if (!pFitness)
			return "Fitness is not of expected type";

		for (float &x : pFitness->getValues())
			x = pGenome->getValues()[0];

		fitness.setCalculated();
		return true;
	}
};

class SingleThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	SingleThreadedPopulationFitnessCalculation(shared_ptr<GenomeFitnessCalculation> genomeFitCalc)
		: m_genomeFitnessCalculation(genomeFitCalc)
	{ 
	}

	~SingleThreadedPopulationFitnessCalculation() { }
	bool_t calculatePopulationFitness(const vector<shared_ptr<Population>> &populations) override;
private:
	shared_ptr<GenomeFitnessCalculation> m_genomeFitnessCalculation;
};

bool_t SingleThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (!m_genomeFitnessCalculation.get())
		return "No genome fitness calculation has been set";

	// First, initialize the calculations
	for (auto &pop : populations)
		for (auto &i : pop->m_individuals)
			if (!i->m_fitness->isCalculated())
				m_genomeFitnessCalculation->startNewCalculation(*i->m_genome);
	
	// Calculate until all is done
	bool allCalculated;
	do
	{
		allCalculated = true;
		for (auto &pop : populations)
		{
			for (auto &i : pop->m_individuals)
			{
				Fitness &f = *i->m_fitness;
				if (!f.isCalculated())
				{
					m_genomeFitnessCalculation->pollCalculate(*i->m_genome, f);
					if (!f.isCalculated())
						allCalculated = false;
				}
			}
		}
	} while (!allCalculated);

	return true;
}

// TODO: set communicator
class MPIPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MPIPopulationFitnessCalculation();
	~MPIPopulationFitnessCalculation();

	// Layout will be exchanged between master and helpers, should be called
	// at same time at master and helpers
	bool_t init(const Genome &referenceGenome,
			    const Fitness &referenceFitness,
				shared_ptr<PopulationFitnessCalculation> &popCalc); // Need reference genome/fitness to be able to serialize it

	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	bool_t calculatePopulationFitness(const vector<shared_ptr<Population>> &populations) override;
	bool_t calculatePopulationFitness_MPIHelper(); // should be called at same time as calculatePopulationFitness
private:
	shared_ptr<Genome> m_referenceGenome;
	shared_ptr<Fitness> m_referenceFitness;
	shared_ptr<Population> m_localPop;
	shared_ptr<PopulationFitnessCalculation> m_localPopulationFitnessCalculation;
};

MPIPopulationFitnessCalculation::MPIPopulationFitnessCalculation()
{
}

MPIPopulationFitnessCalculation::~MPIPopulationFitnessCalculation()
{
}

bool_t MPIPopulationFitnessCalculation::init(const Genome &referenceGenome, const Fitness &referenceFitness,
											 shared_ptr<PopulationFitnessCalculation> &popCalc)
{
	m_referenceGenome = referenceGenome.createCopy(false);
	m_referenceFitness = referenceFitness.createCopy(false);
	m_localPopulationFitnessCalculation = popCalc;
	m_localPop = make_shared<Population>();

	bool_t r = m_referenceGenome->MPI_BroadcastLayout(0, MPI_COMM_WORLD);
	if (!r)
		return "Error broadcasting genome layout: " + r.getErrorString();
	if (!(r = m_referenceFitness->MPI_BroadcastLayout(0, MPI_COMM_WORLD)))
		return "Error broadcasting fitness layout: " + r.getErrorString();

	return true;
}

// TODO: send/check genome/fitness type, to make sure helper is using correct types
// TODO: probably need a fitness calculator instance here?
// TODO: check edge case if no genomes need to be calculated
//       helpers should be signalled to end their helper routine
bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";
	if (!m_localPopulationFitnessCalculation.get())
		return "Local fitness calculation not set";
	
	int mpiSize = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	// TODO: check reference fitness type/layout against population?

	m_localPop->m_individuals.clear();
	vector<vector<pair<Genome *, Fitness *>>> helperGenomes(mpiSize); // TODO: use a reusable array
	int nextHelper = 0;

	// Split the work over the helpers
	for (auto pPop : populations)
	{		
		for (auto &i : pPop->m_individuals)
		{
			// TODO: check for first that type is same as reference genome? Or only for debugging?
			// cerr << typeid(*(i->m_genome.get())).name() << " " << typeid(*(m_referenceGenome.get())).name() << endl;

			if (!i->m_fitness->isCalculated())
			{
				if (nextHelper == 0)
					m_localPop->m_individuals.push_back(i);
				else
					helperGenomes[nextHelper].push_back({i->m_genome.get(), i->m_fitness.get()});

				nextHelper = (nextHelper + 1)%mpiSize;
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
	for (int helper = 1 ; helper < mpiSize ; helper++)
	{
		auto &thisHelperGenomes = helperGenomes[helper];
		int numGenomes = thisHelperGenomes.size();
		MPI_Send(&numGenomes, 1, MPI_INT, helper, 0, MPI_COMM_WORLD);

		for (int individual = 0 ; individual < numGenomes ; individual++)
		{
			Genome *pGenome = thisHelperGenomes[individual].first;
			MPI_Request *pReq = getNextRequest();
			auto r = pGenome->MPI_ISend(helper, individual, MPI_COMM_WORLD, pReq);

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
	for (int helper = 0 ; helper < (int)helperGenomes.size() ; helper++)
	{
		auto &partHelperGenoms = helperGenomes[helper];
		for (int individual = 0 ; individual < (int)partHelperGenoms.size() ; individual++)
		{
			Fitness *pFitness = partHelperGenoms[individual].second;
			assert(pFitness);

			pFitness->MPI_IRecv(helper, individual, MPI_COMM_WORLD, getNextRequest());
			pFitness->setCalculated();
			cerr << "Starting receive for " << helper << "," << individual << endl;
		}
	}

	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

	return true;
}

bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness_MPIHelper()
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";
	if (!m_localPopulationFitnessCalculation.get())
		return "Local fitness calculation not set";

	int numGenomes = 0;
	MPI_Recv(&numGenomes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	cerr << "Calculating " << numGenomes << " in helper" << endl;

	// TODO: recycle previously created instances
	vector<pair<Genome *, Fitness *>> helperGenomes;

	m_localPop->m_individuals.clear();
	for (int i = 0 ; i < numGenomes ; i++) 
	{
		auto genome = m_referenceGenome->createCopy(false);
		auto fitness = m_referenceFitness->createCopy(false);
		fitness->setCalculated(false); // Make sure it's not set accidentally, so that we actually calculate it

		m_localPop->m_individuals.push_back(make_shared<Individual>(genome, fitness));
		helperGenomes.push_back({genome.get(), fitness.get()});
	}

	vector<MPI_Request> requests(numGenomes);
	for (int i = 0 ; i < numGenomes ; i++)
	{
		Genome *pGenome = helperGenomes[i].first;
		auto r = pGenome->MPI_IRecv(0, i, MPI_COMM_WORLD, &requests[i]); 
		if (!r)
			return "Error receiving genome in helper: " + r.getErrorString();
	}

	// Wait for everything to be received
	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	cerr << "Received genomes in helper" << endl;
	
	// TODO: make sure that fitness calculation flags are cleared!!
	auto r = m_localPopulationFitnessCalculation->calculatePopulationFitness({ m_localPop });
	if (!r)
		return "Error in local fitness calculation: " + r.getErrorString();

	// Send fitness results back
	for (int individual = 0 ; individual < (int)helperGenomes.size() ; individual++)
	{
		Fitness *pFitness = helperGenomes[individual].second;
		if (!(r = pFitness->MPI_ISend(0, individual, MPI_COMM_WORLD, &requests[individual])))
			return "Error sending back fitness: " + r.getErrorString();
		
		cerr << "Sending back fitness for " << individual << ": " << pFitness->toString() << endl;

		pFitness->setCalculated(false); // Already clear flag again
	}

	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	
	return true;
}

int main_master(int argc, char *argv[])
{
	cerr << "Master" << endl;

	MPIPopulationFitnessCalculation calc;
	shared_ptr<DummyGenomeFitnessCalculation> genomeCalc = make_shared<DummyGenomeFitnessCalculation>();
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>(genomeCalc);

	int numFloatGenome = 16;
	int numFloatFitness = 2;
	FloatVectorGenome genome(numFloatGenome);
	FloatVectorFitness fitness(numFloatFitness);

	auto r = calc.init(genome, fitness, localCalc);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	shared_ptr<Population> pop = make_shared<Population>();
	for (int i = 0 ; i < 32 ; i++)
		pop->m_individuals.push_back(make_shared<Individual>(
					make_shared<FloatVectorGenome>(numFloatGenome, i),
					make_shared<FloatVectorFitness>(numFloatFitness))); 

	r = calc.calculatePopulationFitness({ pop });
	if (!r)
	{
		cerr << "Error in calculatePopulationFitness: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	for (auto &i : pop->m_individuals)
		cout << i->m_genome->toString() << ": " << i->m_fitness->toString() << endl;

	cerr << "Master done\n";
	return 0;
}

int main_helper(int argc, char *argv[])
{
	cerr << "Helper" << endl;

	MPIPopulationFitnessCalculation calc;
	shared_ptr<DummyGenomeFitnessCalculation> genomeCalc = make_shared<DummyGenomeFitnessCalculation>();
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>(genomeCalc);
	
	FloatVectorGenome genome; // should not now exact layout here, just type
	FloatVectorFitness fitness; // same

	auto r = calc.init(genome, fitness, localCalc);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation in helper: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	if (!(r = calc.calculatePopulationFitness_MPIHelper()))
	{
		cerr << "Error running calculatePopulationFitness_MPIHelper: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	cerr << "Helper done\n";
	return 0;
}

int main2(int argc, char *argv[])
{
	int rank, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	printf("%d/%d\n", rank, procs);

	if (rank == 0)
		return main_master(argc, argv);
	return main_helper(argc, argv);
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int r = main2(argc, argv);
	MPI_Finalize();
	return r;
}


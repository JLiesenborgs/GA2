#include <mpi.h>
#include <errut/booltype.h>
#include <memory>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <cstring>
#include <cassert>
#include <cstdio>

// TODO: we're assuming that rank 0 is the master

// TODO: make communicator configurable

// TODO: for genome and fitness, a copyInto function? So that memory doesn't need to be
//       allocated/reallocated often

// TODO: MPI_Comm seems to be a pointer in OpenMPI

using namespace std;
using namespace errut;

class Genome
{
public:
	Genome() { }
	virtual ~Genome() { }
	virtual shared_ptr<Genome> createCopy(bool copyContents = true) const { return nullptr; }
	virtual string toString() const { return "?"; };

	virtual bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) { return "Not implemented"; }

	virtual bool_t MPI_ISend(int dest, int tag, MPI_Comm communicator, MPI_Request *pRequest) const { return "Not implemented"; }
	virtual bool_t MPI_IRecv(int src, int tag, MPI_Comm communicator, MPI_Request *pRequest) { return "Not implemented"; }
};

class FloatValuesGenome : public Genome
{
public:
	FloatValuesGenome(size_t n = 0) : m_values(n) { }
	FloatValuesGenome(size_t n, float initValue) : m_values(n, initValue) { }
	~FloatValuesGenome() { }

	string toString() const override
	{
		stringstream ss;

		ss << "[";
		for (auto x : m_values)
			ss << " " << x;
		ss << " ]";
		return ss.str();
	}

	shared_ptr<Genome> createCopy(bool copyContents = true) const override
	{
		shared_ptr<FloatValuesGenome> g = make_shared<FloatValuesGenome>(m_values.size());
		if (copyContents && m_values.size() > 0)
		{
			assert(m_values.size() == g->m_values.size());
			memcpy(g->m_values.data(), m_values.data(), sizeof(float)*m_values.size());
		}

		return g;
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
private:
	vector<float> m_values;
};

// TODO: enum/flag calculated
class Fitness
{
public:
	Fitness() : m_calculated(false) { }
	virtual ~Fitness() { }
	virtual shared_ptr<Fitness> createCopy(bool copyContents = true) const { return nullptr; }
	virtual string toString() const { return "?"; };
	bool isCalculated() const { return m_calculated; }
	void setCalculated(bool v = true) { m_calculated = v; }
protected:
	bool m_calculated;
};

class FloatValuesFitness : public Fitness
{
public:
	FloatValuesFitness(size_t n = 0) : m_values(n, 0) { }
	~FloatValuesFitness() { }

	// TODO: how to avoid this copy-paste code?
	shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		shared_ptr<FloatValuesFitness> g = make_shared<FloatValuesFitness>(m_values.size());
		if (copyContents && m_values.size() > 0)
		{
			assert(m_values.size() == g->m_values.size());
			memcpy(g->m_values.data(), m_values.data(), sizeof(float)*m_values.size());
		}

		if (copyContents && isCalculated())
			g->setCalculated();

		return g;
	}

	// TODO: same
	string toString() const override
	{
		stringstream ss;

		ss << "[";
		for (auto x : m_values)
			ss << " " << x;
		ss << " ]";
		return ss.str();
	}

private:
	vector<float> m_values;
};

class GenomeFitnessCalculation
{
public:
	GenomeFitnessCalculation() { }
	virtual ~GenomeFitnessCalculation() { }

	// These are to allow a more async version, but by default the sync version is called
	virtual bool_t startNewCalculation(const Genome &genome) { }
	
	virtual bool_t calculate(const Genome &genome, Fitness &fitness) { return "Not implemented in base class"; }
};

class Individual
{
public:
	Individual(shared_ptr<Genome> genome, shared_ptr<Fitness> fitness)
		: m_genome(genome), m_fitness(fitness) { }
//private:
	shared_ptr<Genome> m_genome;
	shared_ptr<Fitness> m_fitness;
};

// Do we need this? Just a typedef perhaps?
class Population
{
public:
	Population() { }
	~Population() { }

	vector<shared_ptr<Individual>> m_individuals;
};

class PopulationFitnessCalculation
{
public:
	PopulationFitnessCalculation() { }
	virtual ~PopulationFitnessCalculation() { }

	virtual bool_t calculatePopulationFitness(const vector<shared_ptr<Population>> &populations) { return "Not implemented in base class"; }
};

class SingleThreadedPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	SingleThreadedPopulationFitnessCalculation() { }
	~SingleThreadedPopulationFitnessCalculation() { }
	bool_t calculatePopulationFitness(const vector<shared_ptr<Population>> &populations) override;
};

bool_t SingleThreadedPopulationFitnessCalculation::calculatePopulationFitness(const vector<shared_ptr<Population>> &populations)
{
	// TODO
	return true;
}

// TODO: set comminicator
class MPIPopulationFitnessCalculation : public PopulationFitnessCalculation
{
public:
	MPIPopulationFitnessCalculation();
	~MPIPopulationFitnessCalculation();

	// Layout will be exchanged between master and helpers
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
	
	int mpiSize = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

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

	// Receive the calculations

	return true;
}

bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness_MPIHelper()
{
	if (!m_referenceGenome.get() || !m_referenceFitness.get())
		return "Reference genome or fitness not set";

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
	return true;
}

int main_master(int argc, char *argv[])
{
	cerr << "Master" << endl;

	MPIPopulationFitnessCalculation calc;
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>();

	int numFloatGenome = 16;
	int numFloatFitness = 2;
	FloatValuesGenome genome(numFloatGenome);
	FloatValuesFitness fitness(numFloatFitness);

	auto r = calc.init(genome, fitness, localCalc);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	shared_ptr<Population> pop = make_shared<Population>();
	for (int i = 0 ; i < 32 ; i++)
		pop->m_individuals.push_back(make_shared<Individual>(
					make_shared<FloatValuesGenome>(numFloatGenome, i),
					make_shared<FloatValuesFitness>(numFloatFitness))); 

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
	shared_ptr<PopulationFitnessCalculation> localCalc = make_shared<SingleThreadedPopulationFitnessCalculation>();

	FloatValuesGenome genome; // should not now exact layout here, just type
	FloatValuesFitness fitness; // same

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


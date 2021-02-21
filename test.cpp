#include <mpi.h>
#include <errut/booltype.h>
#include <memory>
#include <vector>
#include <iostream>
#include <typeinfo>
#include <cstring>
#include <cassert>
#include <cstdio>

using namespace std;
using namespace errut;

class Genome
{
public:
	Genome() { }
	virtual ~Genome() { }
	virtual shared_ptr<Genome> createCopy(bool copyContents = true) const { return nullptr; }
};

class FloatValuesGenome : public Genome
{
public:
	FloatValuesGenome(size_t n = 0) : m_values(n) { }
	FloatValuesGenome(size_t n, float initValue) : m_values(n, initValue) { }
	~FloatValuesGenome() { }

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
private:
	vector<float> m_values;
};

// TODO: enum/flag calculated
class Fitness
{
public:
	Fitness() { }
	virtual ~Fitness() { }
	virtual shared_ptr<Fitness> createCopy(bool copyContents = true) const { return nullptr; }
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

		return g;
	}
private:
	vector<float> m_values;
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

// TODO: set comminicator
class MPIPopulationFitnessCalculation
{
public:
	MPIPopulationFitnessCalculation();
	~MPIPopulationFitnessCalculation();

	bool_t init(const Genome &referenceGenome,
			    const Fitness &referenceFitness); // Need reference genome/fitness to be able to serialize it

	// TODO: all populations should have exactly the same genomes! (ie same number of floats)
	bool_t calculatePopulationFitness(const vector<Population *> &populations);
	bool_t calculatePopulationFitness_MPIHelper(); // should be called at same time as calculatePopulationFitness
private:
	shared_ptr<Genome> m_referenceGenome;
	shared_ptr<Fitness> m_referenceFitness;
};

MPIPopulationFitnessCalculation::MPIPopulationFitnessCalculation()
{
}

MPIPopulationFitnessCalculation::~MPIPopulationFitnessCalculation()
{
}

bool_t MPIPopulationFitnessCalculation::init(const Genome &referenceGenome, const Fitness &referenceFitness)
{
	m_referenceGenome = referenceGenome.createCopy(false);
	m_referenceFitness = referenceFitness.createCopy(false);
	return true;
}

// TODO: send/check genome/fitness type, to make sure helper is using correct types
bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness(const vector<Population *> &populations)
{
	for (auto pPop : populations)
	{
		for (auto &i : pPop->m_individuals)
		{
			// TODO: check for first that type is same as reference genome
			// cerr << typeid(*(i->m_genome.get())).name() << " " << typeid(*(m_referenceGenome.get())).name() << endl;
		}
	}
	return true;
}

bool_t MPIPopulationFitnessCalculation::calculatePopulationFitness_MPIHelper()
{
	return true;
}

int main_master(int argc, char *argv[])
{
	MPIPopulationFitnessCalculation calc;

	int numFloatGenome = 16;
	int numFloatFitness = 2;
	FloatValuesGenome genome(numFloatGenome);
	FloatValuesFitness fitness(numFloatFitness);

	auto r = calc.init(genome, fitness);
	if (!r)
	{
		cerr << "Can't init MPIPopulationFitnessCalculation: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	Population pop;
	for (int i = 0 ; i < 32 ; i++)
		pop.m_individuals.push_back(make_shared<Individual>(
					make_shared<FloatValuesGenome>(numFloatGenome, i),
					make_shared<FloatValuesFitness>(numFloatFitness))); 

	r = calc.calculatePopulationFitness({ &pop });
	if (!r)
	{
		cerr << "Error in calculatePopulationFitness: " << r.getErrorString() << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	cerr << "Master done\n";
	return 0;
}

int main_helper(int argc, char *argv[])
{
	MPIPopulationFitnessCalculation calc;

	FloatValuesGenome genome; // should not now exact layout here, just type
	FloatValuesFitness fitness; // same

	auto r = calc.init(genome, fitness);
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


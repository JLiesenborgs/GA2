#include "nsga2evolver.h"
#include "evolutionaryalgorithm.h"
#include "vectorgenomefitness.h"
#include "calculation.h"
#include "uniformvectorgenomecrossover.h"
#include "vectorgenomeuniformmutation.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "stopcriterion.h"
#include "mersennerandomnumbergenerator.h"
#include <stdexcept>
#include <limits>
#include <cmath>
#include <random>
#include <iostream>

using namespace eatk;
using namespace std;
using namespace errut;

class BaseCalculation : public GenomeFitnessCalculation
{
public:
	BaseCalculation(size_t N, size_t M, size_t popSize, size_t numGen) : m_N(N), m_M(M), m_popSize(popSize), m_numGen(numGen) { }
	size_t getDimension() const { return m_N; }
	size_t getObjectives() const { return m_M; }

	size_t getNumberOfGenerations() const { return m_numGen; }
	size_t getPopulationSize() const { return m_popSize; }

	void getBounds(size_t dim, double &xMin, double &xMax)
	{
		if (dim >= m_N)
			throw runtime_error("Invalid dimension");
		xMin = numeric_limits<double>::quiet_NaN();
		xMax = numeric_limits<double>::quiet_NaN();
		bounds(dim, xMin, xMax);
		if (isnan(xMin) || isnan(xMax))
			throw runtime_error("Bounds not filled in");
	}

	errut::bool_t calculate(const Genome &genome0, Fitness &fitness0) override
	{
		assert(dynamic_cast<const DoubleVectorGenome *>(&genome0));
		assert(dynamic_cast<DoubleVectorFitness *>(&fitness0));
		const DoubleVectorGenome &genome = static_cast<const DoubleVectorGenome &>(genome0);
		DoubleVectorFitness &fitness = static_cast<DoubleVectorFitness &>(fitness0);
		if (genome.getValues().size() != m_N)
			throw runtime_error("Expecting genome of size " + to_string(m_N) + " but got " + to_string(genome.getValues().size()));

		// Genome: values from -1 to 1, translate to other range
		vector<double> v = genome.getValues();
		for (size_t idx = 0 ; idx < m_N ; idx++)
		{
			double xMin, xMax;
			getBounds(idx, xMin, xMax);
			double frac = (v[idx] + 1.0)/2.0;
			v[idx] = frac*(xMax-xMin) + xMin;
		}

		auto fv = calculate(genome.getValues());
		fitness.setValues(fv);

		if (fitness.getValues().size() != m_M)
			throw runtime_error("Expecting " + to_string(m_M) + " objective values but got " + to_string(fitness.getValues().size()));

		return true;
	}

protected:
	virtual vector<double> calculate(const vector<double> &x) = 0;
	virtual void bounds(size_t dim, double &xMin, double &xMax) = 0;
private:
	size_t m_N, m_M, m_popSize, m_numGen;
};

class SCH : public BaseCalculation
{
public:
	SCH() : BaseCalculation(1, 2, 64, 1000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double x0 = x[0];
		double f1 = x0 * x0;
		double f2 = (x0 - 2.0) * (x0 - 2.0);
		return {f1, f2};
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -1000.0;
		xMax = 1000.0;
	}
};

// TODO: copied this for now, put in common file
class MyEA : public EvolutionaryAlgorithm
{
public:
	MyEA(bool alwaysRecalc = false) : m_alwaysRecalc(alwaysRecalc)
	{
		m_dumpPop = (getenv("DUMPPOP")) ? true : false;
	}

private:
	bool_t onBeforeFitnessCalculation(size_t generation, const std::shared_ptr<Population> &population) override
	{
		if (m_alwaysRecalc)
		{
			for (auto &i : population->individuals())
				i->fitness()->setCalculated(false);
		}
		return true;
	}

	bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<Population> &population) override
	{
		if (m_dumpPop)
		{
			cerr << "Generation " << generation << ":" << endl;
			for (size_t i = 0; i < population->individuals().size(); i++)
				cerr << "  [" << i << "] = " << population->individual(i)->toString() << endl;
			cerr << endl;
		}
		return true;
	}

	bool m_dumpPop;
	bool m_alwaysRecalc;
};

class IndCreation : public IndividualCreation
{
public:
	IndCreation(BaseCalculation &calc, const shared_ptr<RandomNumberGenerator> &rng)
		: m_rng(rng)
	{
		m_N = calc.getDimension();
		m_M = calc.getObjectives();
	}

	shared_ptr<Genome> createInitializedGenome() override
	{
		auto g = make_shared<DoubleVectorGenome>(m_N);
		for (size_t i = 0 ; i < m_N ; i++)
			g->setValue(m_rng->getRandomDouble(-1.0,1.0), i);
		return g;
	}
	
	shared_ptr<Fitness> createEmptyFitness() override
	{
		auto f = make_shared<DoubleVectorFitness>(m_M);
		cout << "f hasRealValues " << f->hasRealValues() << endl;
		return f;
	}

	size_t m_N, m_M;
	shared_ptr<RandomNumberGenerator> m_rng;
};

int main(void)
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));
	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);

	auto problem = make_shared<SCH>();
	MyEA ea;

	double mutFrac = 0.5/(double)problem->getDimension();
	IndCreation creation(*problem, rng);
	NSGA2Evolver evolver(rng,
		make_shared<UniformVectorGenomeCrossover<double>>(rng, false),
		make_shared<VectorGenomeUniformMutation<double>>(mutFrac, -1.0, 1.0, rng),
		make_shared<VectorFitnessComparison<double>>(),
		problem->getObjectives()
		);
	SingleThreadedPopulationFitnessCalculation calc(problem);
	FixedGenerationsStopCriterion stop(problem->getNumberOfGenerations());
	size_t popSize = problem->getPopulationSize();
	
	bool_t r;
	cout << "popSize = " << popSize << endl;
	if (!(r = ea.run(creation, evolver, calc, stop, popSize, popSize, popSize*2)))
		throw runtime_error("Error running EA: " + r.getErrorString());

	auto best = evolver.getBestIndividuals();
	cout << "# Best: " << endl;
	for (auto i : best)
		cout << i->toString() << endl;

	return 0;
}

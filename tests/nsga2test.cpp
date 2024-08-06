#include "nsga2evolver.h"
#include "evolutionaryalgorithm.h"
#include "vectorgenomefitness.h"
#include "calculation.h"
#include "vectorgenomeuniformmutation.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "stopcriterion.h"
#include "mersennerandomnumbergenerator.h"
#include "vectorgenomedelikecrossover.h"
#include "vectorgenomesimulatedbinarycrossover.h"
#include "testfunctions.h"
#include <stdexcept>
#include <limits>
#include <cmath>
#include <random>
#include <iostream>
#include <map>

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

	std::vector<double> getVariableVectorFromGenome(const Genome &genome0)
	{
		assert(dynamic_cast<const DoubleVectorGenome *>(&genome0));
		const DoubleVectorGenome &genome = static_cast<const DoubleVectorGenome &>(genome0);

		if (genome.getValues().size() != m_N)
			throw runtime_error("Expecting genome of size " + to_string(m_N) + " but got " + to_string(genome.getValues().size()));

		return genome.getValues();
	}

	errut::bool_t calculate(const Genome &genome0, Fitness &fitness0) override
	{
		assert(dynamic_cast<DoubleVectorFitness *>(&fitness0));
		DoubleVectorFitness &fitness = static_cast<DoubleVectorFitness &>(fitness0);
		
		auto v = getVariableVectorFromGenome(genome0);
		double penalty = 0;

		for (size_t i = 0 ; i < v.size() ; i++)
		{
			double xMin, xMax;
			double diff = 0;
			getBounds(i, xMin, xMax);
			if (v[i] < xMin)
				diff = xMin-v[i];
			else if (v[i] > xMax)
				diff = v[i] - xMax;

			if (diff)
				penalty += 1e6 + diff;
		}

		if (penalty == 0)
		{
			auto fv = calculate(v);
			fitness.getValues() = fv;
		}
		else
		{
			vector<double> fv(m_M, penalty);
			fitness.getValues() = fv;
		}

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

template <class T>
class BaseCalculationTemplate : public BaseCalculation
{
public:
	template <typename ... Args>
	BaseCalculationTemplate(size_t dim, size_t numObj, size_t popSize, size_t numGen,
							Args&& ... args)
		: BaseCalculation(dim, numObj, popSize, numGen), m_f(std::forward<Args>(args)...) { }
	
	vector<double> calculate(const vector<double> &x) override { return m_f.calculate(x); }
	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		auto [lower,upper] = m_f.getInitialParameterRange();
		xMin = lower.at(dim);
		xMax = upper.at(dim);
	}
private:
	T m_f;
};

class Rosenbrock : public BaseCalculationTemplate<eatk::testfunctions::Rosenbrock>
{
public:
	Rosenbrock() : BaseCalculationTemplate(2, 1, 64, 2000, pair(-2.048, 2.048)) { }
};

class Griewank : public BaseCalculationTemplate<eatk::testfunctions::Griewank>
{
public:
	Griewank() : BaseCalculationTemplate(10, 1, 64, 2000, 10, pair(-400.0, 400.0)) { }
};

class SCH : public BaseCalculationTemplate<eatk::testfunctions::Schaffer>
{
public:
	SCH() : BaseCalculationTemplate(1, 2, 64, 2000, pair(-1000.0, 1000.0)) { }
};

class FON : public BaseCalculationTemplate<eatk::testfunctions::FonsecaFleming>
{
public:
	FON() : BaseCalculationTemplate(3, 2, 64, 2000, pair(-4.0, 4.0)) { }
};

class POL : public BaseCalculationTemplate<eatk::testfunctions::Poloni>
{
public:
	POL() : BaseCalculationTemplate(2, 2, 64, 2000, pair(-M_PI, M_PI)) { }
};

class KUR : public BaseCalculationTemplate<eatk::testfunctions::Kursawe>
{
public:
	KUR() : BaseCalculationTemplate(3, 2, 64, 2000, pair(-5.0, 5.0)) { }
};

class ZDT1 : public BaseCalculationTemplate<eatk::testfunctions::ZitzlerDebThiele1>
{
public:
	ZDT1() : BaseCalculationTemplate(30, 2, 128, 2000, pair(0.0, 1.0)) { }
};

class ZDT2 : public BaseCalculationTemplate<eatk::testfunctions::ZitzlerDebThiele2>
{
public:
	ZDT2() : BaseCalculationTemplate(30, 2, 256, 2000, pair(0.0, 1.0)) { }
};

class ZDT3 : public BaseCalculationTemplate<eatk::testfunctions::ZitzlerDebThiele3>
{
public:
	ZDT3() : BaseCalculationTemplate(30, 2, 128, 2000, pair(0.0, 1.0)) { }
};

class ZDT4 : public BaseCalculationTemplate<eatk::testfunctions::ZitzlerDebThiele4>
{
public:
	ZDT4() : BaseCalculationTemplate(10, 2, 128, 2000, 
		vector{0.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0 },
		vector{1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 }) { }
};

class ZDT6 : public BaseCalculationTemplate<eatk::testfunctions::ZitzlerDebThiele6>
{
public:
	ZDT6() : BaseCalculationTemplate(10, 2, 128, 2000, pair(0.0, 1.0)) { }
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
	IndCreation(const std::shared_ptr<BaseCalculation> &calc, const shared_ptr<RandomNumberGenerator> &rng)
		: m_problem(calc), m_rng(rng)
	{
		m_N = calc->getDimension();
		m_M = calc->getObjectives();
	}

	shared_ptr<Genome> createInitializedGenome() override
	{
		auto g = make_shared<DoubleVectorGenome>(m_N);
		for (size_t i = 0 ; i < m_N ; i++)
		{
			double xMin, xMax;
			m_problem->getBounds(i, xMin, xMax);
			g->getValues()[i] = m_rng->getRandomDouble(xMin, xMax);
		}
		return g;
	}
	
	shared_ptr<Fitness> createEmptyFitness() override
	{
		auto f = make_shared<DoubleVectorFitness>(m_M);
		cout << "f hasRealValues " << f->hasRealValues() << endl;
		return f;
	}

	size_t m_N, m_M;
	shared_ptr<BaseCalculation> m_problem;
	shared_ptr<RandomNumberGenerator> m_rng;
};

std::string getIndividualString(const shared_ptr<Individual> &ind, const shared_ptr<BaseCalculation> &problem)
{
	string s = "genome: ";

	for (auto v : problem->getVariableVectorFromGenome(ind->genomeRef()))
		s += " " + to_string(v);
	s += " | fitness: " + ind->fitness()->toString();
	return s;
}

class Stop : public FixedGenerationsStopCriterion
{
public:
	Stop(size_t numGen, const shared_ptr<BaseCalculation> &problem)
		: FixedGenerationsStopCriterion(numGen), m_problem(problem) { }
	
	errut::bool_t analyze(const PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop) override
	{
		cerr << "# Generation " << generationNumber << " best:" << endl;
		for (auto &i : evolver.getBestIndividuals())
			cerr << getIndividualString(i, m_problem) << endl;
		cerr << endl;
		return FixedGenerationsStopCriterion::analyze(evolver, generationNumber, shouldStop);
	}

	shared_ptr<BaseCalculation> m_problem;
};

int mainCxx(const vector<string> &args)
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));
	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);

	map<string, shared_ptr<BaseCalculation>> problems = {
		{ "Rosenbrock", make_shared<Rosenbrock>() },
		{ "Griewank", make_shared<Griewank>() },
		{ "SCH", make_shared<SCH>() },
		{ "FON", make_shared<FON>() },
		{ "POL", make_shared<POL>() },
		{ "KUR", make_shared<KUR>() },
		{ "ZDT1", make_shared<ZDT1>() },
		{ "ZDT2", make_shared<ZDT2>() },
		{ "ZDT3", make_shared<ZDT3>() },
		{ "ZDT4", make_shared<ZDT4>() },
		{ "ZDT6", make_shared<ZDT6>() }
	};

	auto usage = [&problems]() {
		cerr << "Specify exactly one problem" << endl;
		cerr << "Valid names are:";
		for (auto kv : problems)
			cerr << " " << kv.first;
		cerr << endl;
		return -1;
	};

	if (args.size() != 2)
		return usage();

	auto it = problems.find(args[1]);
	if (it == problems.end())
		return usage();

	auto problem = it->second;
	MyEA ea;

	bool extraParent = false;
	if (getenv("EXTRAPARENT"))
		extraParent = true;
	
	float F = numeric_limits<float>::quiet_NaN();
	float CR = numeric_limits<float>::quiet_NaN();
	if (getenv("F"))
		F = stof(getenv("F"));
	if (getenv("CR"))
		CR = stof(getenv("CR"));

	cerr << "F = " << F << " CR = " << CR << endl;

	IndCreation creation(problem, rng);
	NSGA2Evolver evolver(rng,
		make_shared<VectorGenomeDELikeCrossOver<double>>(rng, extraParent, F, CR),
		//make_shared<VectorGenomeSimulatedBinaryCrossover<double>>(rng, 2.0),
		nullptr,
		make_shared<VectorFitnessComparison<double>>(), problem->getObjectives());

	SingleThreadedPopulationFitnessCalculation calc(problem);
	Stop stop(problem->getNumberOfGenerations(), problem);
	size_t popSize = problem->getPopulationSize();
	
	bool_t r;
	cout << "popSize = " << popSize << endl;
	if (!(r = ea.run(creation, evolver, calc, stop, popSize, popSize, popSize*2)))
		throw runtime_error("Error running EA: " + r.getErrorString());

/*
	auto best = evolver.getBestIndividuals();
	cout << "# Best: " << endl;
	for (auto i : best)
		cout << getIndividualString(i, problem) << endl;
*/
	return 0;
}

int main(int argc, char *argv[])
{
	vector<string> args;
	for (int i = 0 ; i < argc ; i++)
		args.push_back(argv[i]);
	return mainCxx(args);
}

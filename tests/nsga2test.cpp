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

class Rosenbrock : public BaseCalculation
{
public:
	Rosenbrock() : BaseCalculation(2, 1, 64, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);
		double x1 = x[0];
		double x2 = x[1];

		return { (100.0 * (x1*x1 - x2) * (x1*x1 - x2) + (1.0 - x1)*(1.0 - x1)) };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -2.048;
		xMax = 2.048;
	}
};

class Griewank : public BaseCalculation
{
public:
	Griewank() : BaseCalculation(10, 1, 64, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		assert(x.size() == 10);
		double s = 1.0;
        for (size_t i = 0 ; i < 10 ; i++)
            s+= x[i]*x[i]/4000.0;

        double p = 1.0;
        for (size_t i = 0 ; i < 10 ; i++)
            p *= std::cos(x[i]/std::sqrt(i+1));
        
        return { s - p };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -400.0;
		xMax = 400.0;
	}
};

class SCH : public BaseCalculation
{
public:
	SCH() : BaseCalculation(1, 2, 64, 2000) { }
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

class FON : public BaseCalculation
{
public:
	FON() : BaseCalculation(3, 2, 64, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = 0, f2 = 0;
		for (size_t i = 0 ; i < 3 ; i++)
		{
			f1 += pow((x[i] - 1.0/sqrt(3.0)), 2);
			f2 += pow((x[i] + 1.0/sqrt(3.0)), 2);
		}

		f1 = 1.0-exp(-f1);
		f2 = 1.0-exp(-f2);
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -4.0;
		xMax = 4.0;
	}
};

class POL : public BaseCalculation
{
public:
	POL() : BaseCalculation(2, 2, 64, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double A1 = 0.5*sin(1.0) - 2.0*cos(1.0) + sin(2.0) - 1.5*cos(2.0);
		double A2 = 1.5*sin(1.0) - cos(1.0) + 2.0*sin(2.0) - 0.5*cos(2.0);
		double B1 = 0.5*sin(x[0]) - 2.0*cos(x[0]) + sin(x[1]) - 1.5*cos(x[1]);
		double B2 = 1.5*sin(x[0]) - cos(x[0]) + 2.0*sin(x[1]) - 0.5*cos(x[1]);
		double f1 = 1.0 + pow(A1-B1,2) + pow(A2-B2,2);
		double f2 = pow(x[0]+3.0, 2) + pow(x[1]+1.0, 2);
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -M_PI;
		xMax = M_PI;
	}
};

class KUR : public BaseCalculation
{
public:
	KUR() : BaseCalculation(3, 2, 64, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = 0, f2 = 0;
		for (size_t i = 0 ; i < 2 ; i++)
			f1 += -10.0*exp(-0.2*sqrt(x[i]*x[i] + x[i+1]*x[i+1]));
		for (size_t i = 0 ; i < 3 ; i++)
			f2 += std::pow(std::abs(x[i]), 0.8) + 5.0*sin(x[i]*x[i]*x[i]);
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = -5.0;
		xMax = 5.0;
	}
};

class ZDT1 : public BaseCalculation
{
public:
	ZDT1() : BaseCalculation(30, 2, 128, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = x[0];
		double g = 0;
		for (size_t i = 1 ; i < x.size() ; i++)
			g += x[i];

		g = 1.0 + 9.0*g/(double)(x.size() - 1);
		double f2 = g*(1.0-sqrt(x[0]/g));
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = 0.0;
		xMax = 1.0;
	}
};

class ZDT2 : public BaseCalculation
{
public:
	ZDT2() : BaseCalculation(30, 2, 256, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = x[0];
		double g = 0;
		for (size_t i = 1 ; i < x.size() ; i++)
			g += x[i];

		g = 1.0 + 9.0*g/(double)(x.size() - 1);
		double f2 = g*(1.0-pow(x[0]/g, 2));
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = 0.0;
		xMax = 1.0;
	}
};

class ZDT3 : public BaseCalculation
{
public:
	ZDT3() : BaseCalculation(30, 2, 128, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = x[0];
		double g = 0;
		for (size_t i = 1 ; i < x.size() ; i++)
			g += x[i];

		g = 1.0 + 9.0*g/(double)(x.size() - 1);
		double f2 = g*(1.0 - sqrt(x[0]/g) - (x[0]/g)*sin(10.0*M_PI*x[0]));
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = 0.0;
		xMax = 1.0;
	}
};

class ZDT4 : public BaseCalculation
{
public:
	ZDT4() : BaseCalculation(10, 2, 128, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = x[0];
		double g = 0.0;
		for (size_t i = 1 ; i < x.size() ; i++)
			g += pow(x[i], 2) - 10.0*cos(4.0*M_PI*x[i]);

		g = 1.0 + 10.0*((double)x.size()-1.0) + g;
		double f2 = g*(1.0-sqrt(x[0]/g));
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		if (dim == 0)
		{
			xMin = 0.0;
			xMax = 1.0;
		}
		else
		{
			xMin = -5.0;
			xMax = 5.0;
		}
	}
};

class ZDT6 : public BaseCalculation
{
public:
	ZDT6() : BaseCalculation(10, 2, 128, 2000) { }
protected:
	vector<double> calculate(const vector<double> &x) override
	{
		double f1 = 1.0 - exp(-4.0*x[0])*pow(sin(6.0*M_PI*x[0]), 6);
		double g = 0.0;
		for (size_t i = 1 ; i < x.size() ; i++)
			g += x[i];

		g = 1.0 + 9.0*pow(g/((double)x.size() - 1.0), 0.25);
		double f2 = g*(1.0-pow(f1/g, 2));
		return { f1, f2 };
	}

	void bounds(size_t dim, double &xMin, double &xMax) override
	{
		xMin = 0.0;
		xMax = 1.0;
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
		//make_shared<VectorGenomeSimulatedBinaryCrossover<double>>(rng, 1.0),
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

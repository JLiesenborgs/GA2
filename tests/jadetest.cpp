#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "jadeevolver.h"
#include "singlethreadedpopulationfitnesscalculation.h"
#include "testfunctions.h"
#include <random>

using namespace std;
using namespace eatk;
using namespace errut;

template<class VT, class FT>
class VectorValueFitnessCalculation : public GenomeFitnessCalculation
{
public:
	VectorValueFitnessCalculation() : m_evaluations(0) { };

	errut::bool_t calculate(const Genome &genome, Fitness &fitness)
	{
		const VectorGenome<VT> &vg = static_cast<const VectorGenome<VT> &>(genome);
		ValueFitness<FT> &vf = static_cast<ValueFitness<FT> &>(fitness);

		m_evaluations++;

		vf.setValue(calculate(vg.getValues()));
		return true;
	}

	virtual FT calculate(const vector<VT> &x) = 0;
protected:
	size_t m_evaluations;
};

class BaseCalculation : public VectorValueFitnessCalculation<double,double>
{
public:
	size_t getNumberOfEvaluations() const { return m_evaluations; }
	void resetEvaluationCount() { m_evaluations = 0; }
};

template <class T>
class TestFunctionTemplate : public BaseCalculation
{
public:
	template<typename... Args>
	TestFunctionTemplate(Args&&... args) : m_f(std::forward<Args>(args)...) { }

	double calculate(const vector<double> &x) override { return m_f.calculate1(x); } 

	T m_f;
};

typedef TestFunctionTemplate<testfunctions::Sphere> f1_Sphere;
typedef TestFunctionTemplate<testfunctions::Schwefel_2_22> f2;
typedef TestFunctionTemplate<testfunctions::Schwefel_1_2> f3;
typedef TestFunctionTemplate<testfunctions::Schwefel_2_21> f4;
typedef TestFunctionTemplate<testfunctions::GeneralizedRosenbrock> f5;
typedef TestFunctionTemplate<testfunctions::Step_2_Function> f6;
typedef TestFunctionTemplate<testfunctions::QuarticWithNoise> f7; // Doesn't seem to work well
typedef TestFunctionTemplate<testfunctions::ModifiedSchwefel_2_26> f8;
typedef TestFunctionTemplate<testfunctions::Rastrigin> f9;
typedef TestFunctionTemplate<testfunctions::AckleyFunction1> f10;
typedef TestFunctionTemplate<testfunctions::Griewank> f11;
typedef TestFunctionTemplate<testfunctions::GeneralizedPenalizedFunction1> f12;
typedef TestFunctionTemplate<testfunctions::GeneralizedPenalizedFunction2> f13;
typedef TestFunctionTemplate<testfunctions::Branin> f14_Branin;
typedef TestFunctionTemplate<testfunctions::GoldsteinPrice> f15_GoldsteinPrice;
typedef TestFunctionTemplate<testfunctions::Hartman3D> f16_Hartman3D;
typedef TestFunctionTemplate<testfunctions::Hartman6D> f17_Hartman6D;
typedef TestFunctionTemplate<testfunctions::Shekel> f_Shekel;

template<class T>
class ValueToReachStop : public StopCriterion
{
public:
	ValueToReachStop(T value, size_t maxGenerations) : m_value(value), m_maxGen(maxGenerations),m_numGen(0)
	{
		m_dump = (getenv("DUMPPROGRESS"))?true:false;
	}

	bool_t analyze(const PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop) override
	{
		auto &currentBest = evolver.getBestIndividuals();
		if (currentBest.size() != 1)
			return "Expecting current best size 1";
		
		if (m_dump)
			cerr << generationNumber << " " << currentBest[0]->toString() << endl;

		m_numGen = generationNumber;
		const ValueFitness<T> &vf = static_cast<const ValueFitness<T> &>(currentBest[0]->fitnessRef());
		if (vf.getValue() <= m_value)
		{
			shouldStop = true;
			m_best = currentBest[0]->createCopy();

			if (getenv("DUMPFINAL"))
				cerr << "FOUND: " << m_best->toString() << endl;
		}
		else
		{
			if (generationNumber >= m_maxGen)
			{
				shouldStop = true;
				m_numGen = generationNumber;

				if (getenv("DUMPFINAL"))
					cerr << "NOT FINISHED, CURRENT BEST: " << currentBest[0]->toString() << endl;
			}
		}
		return true;			
	}

	size_t getNumberOfGenerations() const { return m_numGen; }
	bool converged() const { return (m_best.get())?true:false; }
	const Individual &bestSolution() const { return *m_best; }
private:
	const T m_value;
	const size_t m_maxGen;
	size_t m_numGen;
	shared_ptr<Individual> m_best;
	bool m_dump;
};

class AvgStd
{
public:
	AvgStd() : m_sum(0), m_sum2(0), m_count(0) { }
	void process(double x)
	{
		m_sum += x;
		m_sum2 += x*x;
		m_count++;
	}
	string toString() const
	{
		double avg = m_sum/m_count;
		double std = std::sqrt(m_sum2/m_count - avg*avg);
		return to_string(avg) + " (" + to_string(std) + ")";
	}
	size_t getCount() const { return m_count; }
private:
	double m_sum;
	double m_sum2;
	size_t m_count;
};

struct Test
{
	string name;
	size_t popSize;
	vector<double> bottom;
	vector<double> top;
	bool constrained;
	double VTR; // Value to reach
	bool recalc;
	size_t maxGenerations;
	shared_ptr<BaseCalculation> calculator;
};

class MyJADE : public JADEEvolver
{
public:
	MyJADE(const std::shared_ptr<RandomNumberGenerator> &rng,
		const std::shared_ptr<DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<FitnessComparison> &fitComp,
		bool archive)
		: JADEEvolver(rng, mut, cross, fitComp, 0, 0.05, 0.1, archive)
	{
		m_logFCR =  (getenv("LOGFCR"))?true:false;
	}
protected:
	void onMutationCrossoverSettings(double muF, double muCR) const override
	{
		if (m_logFCR)
			cerr << muF << " " << muCR << endl;
	}
private:
	bool m_logFCR;
};

class MyEA : public EvolutionaryAlgorithm
{
public:
	MyEA(bool alwaysRecalc = false) : m_alwaysRecalc(alwaysRecalc)
	{
		m_dumpPop = (getenv("DUMPPOP"))?true:false;
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
			for (size_t i = 0 ; i < population->individuals().size() ; i++)
				cerr << "  [" << i << "] = " << population->individual(i)->toString() << endl;
			cerr << endl;
		}
		return true;
	}

	bool m_dumpPop;
	bool m_alwaysRecalc;
};

int main(int argc, char const *argv[])
{
	string testName;
	if (argc <= 1) { } // nothing to do
	else if (argc == 2)
		testName = argv[1];
	else
	{
		cerr << "Usage: " << argv[0] << " [testname]" << endl;
		return -1;
	}

	bool didRunTest = false;

	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	bool archive = true;
	if (getenv("NOARCH"))
		archive = false;

	cout << "# Seed = " << seed << endl;
	cout << "# Runs = " << numRuns << endl;

	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);

	vector<Test> tests {
		{ "f1_Sphere_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f1_Sphere>(30, pair(-100.0, 100.0)) },
		{ "f1_Sphere_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f1_Sphere>(100, pair(-100.0, 100.0)) },
		{ "f2_30", 100, vector<double>(30, -10), vector<double>(30, 10), false, 1e-8, false, 100000, make_shared<f2>(30, pair(-10.0,10.0)) },
		{ "f2_100", 400, vector<double>(100, -10), vector<double>(100, 10), false, 1e-8, false, 100000, make_shared<f2>(100, pair(-10.0,10.0)) },
		{ "f3_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f3>(30, pair(-100.0,100.0)) },
		{ "f3_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f3>(100, pair(-100.0,100.0)) },
		{ "f4_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f4>(30, pair(-100.0,100.0)) },
		{ "f4_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f4>(100, pair(-100.0,100.0)) },
		{ "f5_30", 100, vector<double>(30, -30), vector<double>(30, 30), false, 1e-8, false, 100000, make_shared<f5>(30, pair(-30.0,30.0)) },
		{ "f5_100", 400, vector<double>(100, -30), vector<double>(100, 30), false, 1e-8, false, 100000, make_shared<f5>(100, pair(-30.0,30.0)) },
		{ "f6_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f6>(30, pair(-100.0,100.0)) },
		{ "f6_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f6>(100, pair(-100.0,100.0)) },
		// Doesn't seem to work??
		// { "f7_30", 100, vector<double>(30, -1.28), vector<double>(30, 1.28), false, 1e-2, true, 100000, make_shared<f7>(30, rng, pair(-1.28,1.28)) },
		// { "f7_100", 400, vector<double>(100, -1.28), vector<double>(100, 1.28), false, 1e-2, true, 100000, make_shared<f7>(100, rng, pair(-1.28,1.28)) },
		{ "f8_30", 100, vector<double>(30, -500), vector<double>(30, 500), true, 1e-8, false, 100000, make_shared<f8>(30, pair(-500.0,500.0)) },
		{ "f8_100", 400, vector<double>(100, -500), vector<double>(100, 500), true, 1e-8, false, 100000, make_shared<f8>(100, pair(-500.0,500.0)) },
		{ "f9_30", 100, vector<double>(30, -5.12), vector<double>(30, 5.12), false, 1e-8, false, 100000, make_shared<f9>(30, pair(-5.12,5.12)) },
		{ "f9_100", 400, vector<double>(100, -5.12), vector<double>(100, 5.12), false, 1e-8, false, 100000, make_shared<f9>(100, pair(-5.12,5.12)) },
		{ "f10_30", 100, vector<double>(30, -32), vector<double>(30, 32), false, 1e-8, false, 100000, make_shared<f10>(30, pair(-32.0,32.0)) },
		{ "f10_100", 400, vector<double>(100, -32), vector<double>(100, 32), false, 1e-8, false, 100000, make_shared<f10>(100, pair(-32.0,32.0)) },
		{ "f11_30", 100, vector<double>(30, -600), vector<double>(30, 600), false, 1e-8, false, 100000, make_shared<f11>(30, pair(-600.0,600.0)) },
		{ "f11_100", 400, vector<double>(100, -600), vector<double>(100, 600), false, 1e-8, false, 100000, make_shared<f11>(100, pair(-600.0,600.0)) },
		{ "f12_30", 100, vector<double>(30, -50), vector<double>(30, 50), false, 1e-8, false, 100000, make_shared<f12>(30, pair(-50.0,50.0)) },
		{ "f12_100", 400, vector<double>(100, -50), vector<double>(100, 50), false, 1e-8, false, 100000, make_shared<f12>(100, pair(-50.0,50.0)) },
		{ "f13_30", 100, vector<double>(30, -50), vector<double>(30, 50), false, 1e-8, false, 100000, make_shared<f13>(30, pair(-50.0,50.0)) },
		{ "f13_100", 400, vector<double>(100, -50), vector<double>(100, 50), false, 1e-8, false, 100000, make_shared<f13>(100, pair(-50.0,50.0)) },
		{ "f14_Branin", 30, { -5.0, 0.0 }, { 10.0, 15.0 }, false,  0.3978873577 + 1e-8, false, 100000, make_shared<f14_Branin>(vector{-5.0,0.0},vector{10.0,15.0}) },
		{ "f15_GoldsteinPrice", 30, { -2.0, -2.0 }, { 2.0, 2.0 }, false,  3.0 + 1e-8, false, 100000, make_shared<f15_GoldsteinPrice>(pair(-2.0,2.0)) },
		{ "f16_Hartman3D", 30, { 0, 0, 0 }, { 1, 1, 1 }, false,  -3.862779787 + 1e-8, false, 100000, make_shared<f16_Hartman3D>(pair(0.0,1.0)) },
		{ "f17_Hartman6D", 30, { 0, 0, 0, 0, 0, 0 }, { 1, 1, 1, 1, 1, 1 }, false,  -3.322368011 + 1e-8, false, 100000, make_shared<f17_Hartman6D>(pair(0.0, 1.0)) },
		{ "f18_Shekel5", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.153199679 + 1e-8, false, 100000, make_shared<f_Shekel>(5, pair(0.0,1.0)) },
		{ "f19_Shekel7", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.402915336 + 1e-8, false, 100000, make_shared<f_Shekel>(7, pair(0.0,1.0)) },
		{ "f20_Shekel10", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.536443153 + 1e-8, false, 100000, make_shared<f_Shekel>(10, pair(0.0,1.0)) },
	};

	auto comp = make_shared<ValueFitnessComparison<double>>();

	for (const auto &test : tests)
	{
		if (testName.size() > 0)
		{
			if (test.name != testName)
				continue;
		}
		didRunTest = true;
		cout << test.name << ": "; 

		AvgStd nfe;
		AvgStd genAvg;
		size_t successCount = 0;
		size_t totalCalc = 0;

		for (size_t run = 0 ; run < numRuns ; run++)
		{
			auto mut = make_shared<VectorDifferentialEvolutionMutation<double>>();
			shared_ptr<VectorDifferentialEvolutionCrossover<double>> cross;
			
			if (test.constrained)
				cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(rng, test.bottom, test.top);
			else
				cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(rng);
			
			VectorDifferentialEvolutionIndividualCreation<double,double> creation(test.bottom, test.top, rng);
			ValueToReachStop<double> stop(test.VTR, test.maxGenerations);
		
			MyJADE evolver(rng, mut, cross, comp, archive);
			SingleThreadedPopulationFitnessCalculation popCalc(test.calculator);

			MyEA ea(test.recalc);
			bool_t r = ea.run(creation, evolver, popCalc, stop, test.popSize, test.popSize, test.popSize*2);
			if (!r)
			{
				cerr << r.getErrorString() << endl;
				return -1;
			}

			if (!stop.converged())
			{
				// cout << "  Failed to converge after " << stop.getNumberOfGenerations() << " generations" << endl;
			}
			else
			{
				//cout << "  Found after " << stop.getNumberOfGenerations() << " generations, " << test.calculator->getNumberOfEvaluations()
				// 	 << " evaluations: " << stop.bestSolution().toString() << endl;

				nfe.process(test.calculator->getNumberOfEvaluations());
				genAvg.process(stop.getNumberOfGenerations());
			}

			test.calculator->resetEvaluationCount();
		}

		cout << nfe.getCount() << "/" << numRuns;
		cout << " , " << genAvg.toString() << " gen";
		cout << " , nfe " << nfe.toString() << " , NP = " << test.popSize << ((archive)?" , archive":", noarchive");
		if (successCount != numRuns)
			cout << " !!";
		cout << endl;
	}

	if (!didRunTest)
	{
		cerr << "No tests executed, available names are: " << endl;
		for (const auto &test : tests)
			cerr << "  " << test.name << endl;
	}

	return 0;
}

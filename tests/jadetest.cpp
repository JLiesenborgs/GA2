#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "jadeevolver.h"
#include "singlethreadedpopulationfitnesscalculation.h"
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

class f1_Sphere : public BaseCalculation
{
public:
	f1_Sphere(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return sumSquared;
	}
private:
	const size_t m_D;
};

class f2 : public BaseCalculation
{
public:
	f2(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0, prod = 1;
		for (auto v : x)
		{
			sum += std::abs(v);
			prod *= std::abs(v);
		}
		return sum + prod;
	}
private:
	const size_t m_D;
};

class f3 : public BaseCalculation
{
public:
	f3(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		
		for (size_t i = 0 ; i < m_D ; i++)
		{
			double sub = 0;
			for (size_t j = 0 ; j <= i ; j++)
				sub += x[j];
			sum += sub*sub;
		}
		return sum;
	}
private:
	const size_t m_D;
};

class f4 : public BaseCalculation
{
public:
	f4(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double max = -numeric_limits<double>::max();
		for (auto v : x)
			if (std::abs(v) > max)
				max = std::abs(v);
		return max;
	}
private:
	const size_t m_D;
};

class f5 : public BaseCalculation
{
public:
	f5(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		for (size_t i = 0 ; i < m_D-1 ; i++)
			sum += 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
		return sum;
	}
private:
	const size_t m_D;
};

class f6 : public BaseCalculation
{
public:
	f6(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		for (auto v : x)
			sum += std::floor(v+0.5)*std::floor(v+0.5);
		return sum;
	}
private:
	const size_t m_D;
};

class f7 : public BaseCalculation // Doesn't seem to work well
{
public:
	f7(size_t D, const shared_ptr<RandomNumberGenerator> &rng) : m_D(D), m_rng(rng) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		for (size_t i = 0 ; i < m_D ; i++)
			sum += (i+1.0)*x[i]*x[i]*x[i]*x[i];
		return sum + m_rng->getRandomDouble();
	}
private:
	const size_t m_D;
	shared_ptr<RandomNumberGenerator> m_rng;
};

class f8 : public BaseCalculation
{
public:
	f8(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);

		double sum = 0;
		for (auto v : x)
			sum += -v*std::sin(std::sqrt(std::abs(v)));
		return sum + 418.98288727243369*m_D;
	}
private:
	const size_t m_D;
};

class f9 : public BaseCalculation
{
public:
	f9(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0;
		for (auto v : x)
			sum += v*v -10.0*std::cos(2.0*M_PI*v) + 10.0;
		return sum;
	}
private:
	const size_t m_D;
};

class f10 : public BaseCalculation
{
public:
	f10(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum1 = 0, sum2 = 0;
		for (auto v : x)
		{
			sum1 += v*v;
			sum2 += std::cos(2.0*M_PI*v);
		}

		return -20.0*std::exp(-0.2*std::sqrt(sum1/m_D)) - std::exp(sum2/m_D) + 20 + std::exp(1);
	}
private:
	const size_t m_D;
};

class f11 : public BaseCalculation
{
public:
	f11(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		double sum = 0, prod = 1;
		for (size_t i = 0 ; i < m_D ; i++)
		{
			sum += x[i]*x[i];
			prod *= std::cos(x[i]/std::sqrt(i+1.0));
		}
		return sum/4000.0 - prod + 1.0;
	}
private:
	const size_t m_D;
};

inline double u(double z, double a, double k, double m)
{
	if (z > a)
		return k*std::pow(z-a, m);
	if (z < -a)
		return k*std::pow(-z-a, m);
	return 0;
};

class f12 : public BaseCalculation
{
public:
	f12(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);
		auto y = [&x](size_t i) { return 1.0+0.25*(x[i] + 1.0); };

		double sum = 10.0*std::sin(M_PI*y(0))*std::sin(M_PI*y(0));

		for (size_t i = 0 ; i < m_D-1 ; i++)
			sum += (y(i)-1.0)*(y(i)-1.0)*(1.0+10.0*std::sin(M_PI*y(i+1))*std::sin(M_PI*y(i+1)));

		sum += (y(m_D-1)-1.0)*(y(m_D-1)-1.0);
		sum *= M_PI/m_D;

		for (size_t i = 0 ; i < m_D ; i++)
			sum += u(x[i], 10, 100, 4);
		return sum;
	}
private:
	const size_t m_D;
};

class f13 : public BaseCalculation
{
public:
	f13(size_t D) : m_D(D) { }
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == m_D);

		double sum = std::sin(3*M_PI*x[0])*std::sin(3*M_PI*x[0]);

		for (size_t i = 0 ; i < m_D-1 ; i++)
			sum += (x[i] - 1.0)*(x[i] - 1.0)*(1.0+std::sin(3*M_PI*x[i+1])*std::sin(3*M_PI*x[i+1]));

		sum += (x[m_D-1]-1.0)*(x[m_D-1]-1.0)*(1.0*std::sin(2*M_PI*x[m_D-1])*std::sin(2*M_PI*x[m_D-1]));
		sum *= 0.1;

		for (size_t i = 0 ; i < m_D ; i++)
			sum += u(x[i], 5, 100, 4);

		return sum;
	}
private:
	const size_t m_D;
};

class f14_Branin : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);
		const double a = 1.0;
		const double b = 5.1/(4.0*M_PI*M_PI);
		const double c = 5.0/M_PI;
		const double r = 6.0;
		const double s = 10.0;
		const double t = 1.0/(8.0*M_PI);

		return a*std::pow(x[1] - b*x[0]*x[0] + c*x[0] - r,2) + s*(1.0-t)*std::cos(x[0]) + s;
	}
};

class f15_GoldsteinPrice : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);

		return (1.0+std::pow(x[0]+x[1]+1.0,2)*(19.0-14.0*x[0]+3.0*x[0]*x[0]-14.0*x[1]+6.0*x[0]*x[1]+3.0*x[1]*x[1]))*(
			30.0+std::pow(2.0*x[0]-3.0*x[1],2)*(18.0-32.0*x[0]+12.0*x[0]*x[0]+48.0*x[1]-36.0*x[0]*x[1]+27.0*x[1]*x[1])
		);
	}
};

class f16_Hartman3D : public BaseCalculation
{
public:
	f16_Hartman3D()
	{
		alpha = { 1.0, 1.2, 3.0, 3.2 };
		A = {
			{ 3.0, 10.0, 30.0},
			{ 0.1, 10.0, 35.0},
			{ 3.0, 10.0, 30.0},
			{ 0.1, 10.0, 35.0}
		};

		P = {
			{ 3689, 1170, 2673 },
			{ 4699, 4387, 7470 },
			{ 1091, 8732, 5547 },
			{  381, 5743, 8828 }
		};
		for (auto &v : P)
			for (auto &x : v)
				x *= 1e-4;
	}

	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 3);
		double sum = 0;
		for (size_t i = 0 ; i < 4 ; i++)
		{
			double s = 0;
			for (size_t j = 0 ; j < 3 ; j++)
				s += A[i][j]*std::pow(x[j]-P[i][j],2);

			sum += alpha[i]*std::exp(-s);
		}
		return -sum;
	}
private:
	vector<double> alpha;
	vector<vector<double>> A, P;
};

class f17_Hartman6D : public BaseCalculation
{
public:
	f17_Hartman6D()
	{
		alpha = { 1.0, 1.2, 3.0, 3.2 };
		A = {
			{ 10, 3, 17, 3.5, 1.7, 8 },
			{ 0.05, 10, 17, 0.1, 8, 14 },
			{ 3, 3.5, 1.7, 10, 17, 8 },
			{ 17, 8, 0.05, 10, 0.1, 14 }
		};

		P = {
			{ 1312, 1696, 5569, 124, 8283, 5886 },
			{ 2329, 4135, 8307, 3736, 1004, 9991 },
			{ 2348, 1451, 3522, 2883, 3047, 6650 },
			{ 4047, 8828, 8732, 5743, 1091, 381 }
		};
		for (auto &v : P)
			for (auto &x : v)
				x *= 1e-4;
	}

	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 6);
		double sum = 0;
		for (size_t i = 0 ; i < 4 ; i++)
		{
			double s = 0;
			for (size_t j = 0 ; j < 6 ; j++)
				s += A[i][j]*std::pow(x[j]-P[i][j],2);

			sum += alpha[i]*std::exp(-s);
		}
		return -sum;
	}
private:
	vector<double> alpha;
	vector<vector<double>> A, P;
};

class f_Shekel : public BaseCalculation
{
public:
	f_Shekel(size_t m) : m_m(m)
	{
		beta = { 1, 2, 2, 4, 4, 6, 3, 7, 5, 5 };
		for (auto &b : beta)
			b *= 0.1;
		
		C = {
			{ 4, 1, 8, 6, 3, 2, 5, 8, 6, 7 },
			{ 4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6 },
			{ 4, 1, 8, 6, 3, 2, 5, 8, 6, 7 },
			{ 4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6 }
		};
	}

	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 4);
		assert(m_m <= 10);
		double sum = 0;
		for (size_t i = 0 ; i < m_m ; i++)
		{
			double s = beta[i];
			for (size_t j = 0 ; j < 4 ; j++)
				s += std::pow(x[j] - C[j][i],2);

			sum += 1.0/s;
		}
		return -sum;
	}
private:
	const size_t m_m;
	vector<double> beta;
	vector<vector<double>> C;
};


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
		{ "f1_Sphere_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f1_Sphere>(30) },
		{ "f1_Sphere_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f1_Sphere>(100) },
		{ "f2_30", 100, vector<double>(30, -10), vector<double>(30, 10), false, 1e-8, false, 100000, make_shared<f2>(30) },
		{ "f2_100", 400, vector<double>(100, -10), vector<double>(100, 10), false, 1e-8, false, 100000, make_shared<f2>(100) },
		{ "f3_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f3>(30) },
		{ "f3_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f3>(100) },
		{ "f4_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f4>(30) },
		{ "f4_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f4>(100) },
		{ "f5_30", 100, vector<double>(30, -30), vector<double>(30, 30), false, 1e-8, false, 100000, make_shared<f5>(30) },
		{ "f5_100", 400, vector<double>(100, -30), vector<double>(100, 30), false, 1e-8, false, 100000, make_shared<f5>(100) },
		{ "f6_30", 100, vector<double>(30, -100), vector<double>(30, 100), false, 1e-8, false, 100000, make_shared<f6>(30) },
		{ "f6_100", 400, vector<double>(100, -100), vector<double>(100, 100), false, 1e-8, false, 100000, make_shared<f6>(100) },
		// Doesn't seem to work??
		// { "f7_30", 100, vector<double>(30, -1.28), vector<double>(30, 1.28), false, 1e-2, true, 100000, make_shared<f7>(30, rng) },
		// { "f7_100", 400, vector<double>(100, -1.28), vector<double>(100, 1.28), false, 1e-2, true, 100000, make_shared<f7>(100, rng) },
		{ "f8_30", 100, vector<double>(30, -500), vector<double>(30, 500), true, 1e-8, false, 100000, make_shared<f8>(30) },
		{ "f8_100", 400, vector<double>(100, -500), vector<double>(100, 500), true, 1e-8, false, 100000, make_shared<f8>(100) },
		{ "f9_30", 100, vector<double>(30, -5.12), vector<double>(30, 5.12), false, 1e-8, false, 100000, make_shared<f9>(30) },
		{ "f9_100", 400, vector<double>(100, -5.12), vector<double>(100, 5.12), false, 1e-8, false, 100000, make_shared<f9>(100) },
		{ "f10_30", 100, vector<double>(30, -32), vector<double>(30, 32), false, 1e-8, false, 100000, make_shared<f10>(30) },
		{ "f10_100", 400, vector<double>(100, -32), vector<double>(100, 32), false, 1e-8, false, 100000, make_shared<f10>(100) },
		{ "f11_30", 100, vector<double>(30, -600), vector<double>(30, 600), false, 1e-8, false, 100000, make_shared<f11>(30) },
		{ "f11_100", 400, vector<double>(100, -600), vector<double>(100, 600), false, 1e-8, false, 100000, make_shared<f11>(100) },
		{ "f12_30", 100, vector<double>(30, -50), vector<double>(30, 50), false, 1e-8, false, 100000, make_shared<f12>(30) },
		{ "f12_100", 400, vector<double>(100, -50), vector<double>(100, 50), false, 1e-8, false, 100000, make_shared<f12>(100) },
		{ "f13_30", 100, vector<double>(30, -50), vector<double>(30, 50), false, 1e-8, false, 100000, make_shared<f13>(30) },
		{ "f13_100", 400, vector<double>(100, -50), vector<double>(100, 50), false, 1e-8, false, 100000, make_shared<f13>(100) },
		{ "f14_Branin", 30, { -5.0, 0.0 }, { 10.0, 15.0 }, false,  0.3978873577 + 1e-8, false, 100000, make_shared<f14_Branin>() },
		{ "f15_GoldsteinPrice", 30, { -2.0, -2.0 }, { 2.0, 2.0 }, false,  3.0 + 1e-8, false, 100000, make_shared<f15_GoldsteinPrice>() },
		{ "f16_Hartman3D", 30, { 0, 0, 0 }, { 1, 1, 1 }, false,  -3.862779787 + 1e-8, false, 100000, make_shared<f16_Hartman3D>() },
		{ "f17_Hartman6D", 30, { 0, 0, 0, 0, 0, 0 }, { 1, 1, 1, 1, 1, 1 }, false,  -3.322368011 + 1e-8, false, 100000, make_shared<f17_Hartman6D>() },
		{ "f18_Shekel5", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.153199679 + 1e-8, false, 100000, make_shared<f_Shekel>(5) },
		{ "f19_Shekel7", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.402915336 + 1e-8, false, 100000, make_shared<f_Shekel>(7) },
		{ "f20_Shekel10", 30, { 0, 0, 0, 0 }, { 10, 10, 10, 10 }, false, -10.536443153 + 1e-8, false, 100000, make_shared<f_Shekel>(10) },
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

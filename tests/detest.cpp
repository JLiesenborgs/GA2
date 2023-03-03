
#include "vectordifferentialevolution.h"
#include "evolutionaryalgorithm.h"
#include "mersennerandomnumbergenerator.h"
#include "valuefitness.h"
#include "differentialevolutionevolver.h"
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
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 3);
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return sumSquared;
	}
};

class f2_Rosenbrock : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);
		double x1 = x[0];
		double x2 = x[1];

		return (100.0 * (x1*x1 - x2) * (x1*x1 - x2) + (1.0 - x1)*(1.0 - x1));
	}
};

class f3_Mod3rdDeJong : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 5);
		double r = 30.0;
		for (auto v : x)
		{
			if (std::abs(v) < 5.12)
				r += std::floor(v);
			else if (v > 5.12)
				r += 30.0*(v - 5.12);
			else
				r += 30.0*(5.12 - v);
		}
		return r;
	}
};

class f4_Mod4thDeJong_quartic : public BaseCalculation
{
public:
	f4_Mod4thDeJong_quartic(const shared_ptr<RandomNumberGenerator> &rng) : m_rng(rng) { }
	static vector<double> lower() { return vector<double>(30, -1.28); }
	static vector<double> upper() { return vector<double>(30, 1.28); }

	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 30);
		double r = 0;
	
		for (size_t i = 0 ; i < 30 ; i++)
		{
			double j = i+1.0;
			double eta = m_rng->getRandomDouble();
			double xSquared = x[i]*x[i];
			double xFourth = xSquared*xSquared;
			r += j*xFourth + eta;
		}
		return r;
	};
private:
	shared_ptr<RandomNumberGenerator> m_rng;
};

class f5_Foxholes : public BaseCalculation
{
public:
	f5_Foxholes()
	{
		a = { 
            { -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32,-32, -16, 0, 16, 32,-32, -16, 0, 16, 32},
            { -32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0,0,0,0,0, 16,16,16,16,16, 32,32,32,32,32 }      
        };
	}

	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);
        double s = 0.002;
        for (size_t i = 0 ; i < 25 ; i++)
        {
            double term = i + 1.0;

            double diff0 = x[0] - a[0][i];
            double diff1 = x[1] - a[1][i];

            term += diff0*diff0*diff0*diff0*diff0*diff0;
            term += diff1*diff1*diff1*diff1*diff1*diff1;

            s += 1.0/term;
        }
        return 1.0/s;
	}
private:
    vector<vector<double>> a;
};

class f6_Corana : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 4);

		static const double d[] = {1,1000,10,100};
        auto sign = [](double x) -> double
        {
            if (x > 0)
                return 1.0;
            if (x < 0)
                return -1.0;
            return 0.0;
        };

        double s = 0;
        for (size_t i = 0 ; i < 4 ; i++)
        {
            double z = std::floor(std::abs(x[i]/0.2) + 0.49999)*sign(x[i])*0.2;
            if (std::abs(x[i] - z) < 0.05)
                s += 0.15*(z-0.05*sign(z))*(z-0.05*sign(z))*d[i];
            else
                s += d[i]*x[i]*x[i];
        }
        return s;
	}
};

class f7_Griewangk : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 10);
		double s = 1.0;
        for (size_t i = 0 ; i < 10 ; i++)
            s+= x[i]*x[i]/4000.0;

        double p = 1.0;
        for (size_t i = 0 ; i < 10 ; i++)
            p *= std::cos(x[i]/std::sqrt(i+1));
        
        return s - p;
	}
};

class f8_Zimmermann : public BaseCalculation
{
public:
	double calculate(const vector<double> &x) override
	{
		assert(x.size() == 2);
        double f = 9.0 - x[0] - x[1];
        double constraint1 = (x[0] - 3.0)*(x[0] - 3.0) + (x[1] - 2.0)*(x[1] - 2.0);
        double constraint2 =  x[0]*x[1];

        if (constraint1 > 16)
            f += 100 + 100*(constraint1 - 16);
        if (constraint2 > 14)
            f += 100 + 100*(constraint2 - 14);
        return f;
	}
};

class f9_k4_Poly : public BaseCalculation
{
private:
    vector<double> z;
public:
    f9_k4_Poly()
    {
        z.push_back(-1.2);
        for (size_t i = 0 ; i < 60 ; i++)
            z.push_back(i*2.0/(60-1) + (-1.0));
        z.push_back(1.2);
    }

	double calculate(const vector<double> &x) override
    {
        auto T8 = [](double z) 
        {
            if (z == 1.2 || z == -1.2)
                return 72.6606669;

            double z2 = z*z;
            double z4 = z2*z2; 
            double z6 = z4*z2;
            double z8 = z4*z4;
            return 1.0 - 32.0*z2 + 160.0*z4 -256.0*z6 + 128.0*z8;
        };

        auto f9 = [](const vector<double> &x, double z)
        {
            double s = 0;
            double zj = 1.0;
            for (auto v : x)
            {
                s += v*zj;
                zj *= z;
            }
            return s;
        };

        double sumDiff = 0;
        for (auto zz : z)
        {
            double pred = f9(x, zz);
            double real = T8(zz);
            double diff = (pred-real);
            double diffSquared = diff*diff;
            sumDiff += diffSquared;
        }
        return sumDiff;
    }
};

class f9_k8_Poly : public BaseCalculation
{
private:
    vector<double> z;
public:
    f9_k8_Poly()
    {
        z.push_back(-1.2);
        for (size_t i = 0 ; i < 60 ; i++)
            z.push_back(i*2.0/(60-1) + (-1.0));
        z.push_back(1.2);
    }

	double calculate(const vector<double> &x) override
    {
        auto T16 = [](double z) 
        {
            if (z == 1.2 || z == -1.2)
                return 10558.1450229;

            double z2 = z*z;
            double z4 = z2*z2; 
            double z6 = z4*z2;
            double z8 = z4*z4;
            double z10 = z4*z6;
            double z12 = z6*z6;
            double z14 = z8*z6;
            double z16 = z8*z8;
            return 1.0 - 128.0*z2 + 2688.0*z4 -21504.0*z6 + 84480.0*z8
                  -180224.0*z10 + 212992.0*z12 -131072.0*z14 + 32768.0*z16;
        };

        auto f9 = [](const vector<double> &x, double z)
        {
            double s = 0;
            double zj = 1.0;
            for (auto v : x)
            {
                s += v*zj;
                zj *= z;
            }
            return s;
        };

        double sumDiff = 0;
        for (auto zz : z)
        {
            double pred = f9(x, zz);
            double real = T16(zz);
            double diff = (pred-real);
            double diffSquared = diff*diff;
            sumDiff += diffSquared;
        }
        return sumDiff;
    }
};

class f11_HyperEllipsoid : public BaseCalculation
{
private:
    const size_t m_D;
public:
    f11_HyperEllipsoid(size_t D) : m_D(D) { }

	double calculate(const vector<double> &x) override
    {
        double s = 0;

        for (size_t j = 0 ; j < m_D ; j++)
            s += (j+1.0)*(j+1.0)*x[j]*x[j];

        return s;
    }
};

class f12_Katsuura : public BaseCalculation
{
private:
    const size_t m_D;
public:
    f12_Katsuura(size_t D) : m_D(D) { }

    double calculate(const vector<double> &x) override
    {
        double p = 1.0;

        for (size_t j = 0 ; j < m_D ; j++)
        {
            double s = 0.0;
            double twok = 2.0;
            for (size_t k = 1 ; k < 33 ; k++)
            {
                s += std::floor(std::abs(twok*x[j]))/twok;
                twok *= 2.0;
            }

            s *= (j+1.0);
            s += 1.0;

            p *= s;
        }
        return p;
    }
};

template<class T>
class ValueToReachStop : public StopCriterion
{
public:
	ValueToReachStop(T value, size_t maxGenerations) : m_value(value), m_maxGen(maxGenerations),m_numGen(0)
	{
		m_dump = (getenv("DUMPPROGRESS"))?true:false;
	}

	bool_t analyze(const std::vector<std::shared_ptr<Individual>> &currentBest, size_t generationNumber, bool &shouldStop)
	{
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

struct Test
{
	string name;
	size_t popSize;
	double F;
	double CR;
	vector<double> bottom;
	vector<double> top;
	double VTR; // Value to reach
	size_t maxGenerations;
	shared_ptr<BaseCalculation> calculator;
};

class MyEA : public EvolutionaryAlgorithm
{
public:
	MyEA()
	{
		m_dumpPop = (getenv("DUMPPOP"))?true:false;
	}
private:
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
};

int main(int argc, char const *argv[])
{
	random_device rndDev;
	unsigned int seed = rndDev();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	cout << "# Seed = " << seed << endl;
	cout << "# Runs = " << numRuns << endl;

	shared_ptr<RandomNumberGenerator> rng = make_shared<MersenneRandomNumberGenerator>(seed);

	vector<Test> tests {
		{ "f1_Sphere", 10, 0.5, 0.3, { -5.12, -5.12, -5.12 }, { 5.12, 5.12, 5.12 }, 1e-6, 100000, make_shared<f1_Sphere>() },
		{ "f2_Rosenbrock", 10, 0.9, 0.9, { -2.048, -2.048 }, { 2.048, 2.048 }, 1e-6, 100000, make_shared<f2_Rosenbrock>() },
		{ "f3_Mod3rdDeJong_step", 10, 0.9, 0, { -5.12, -5.12, -5.12, -5.12, -5.12 }, { 5.12, 5.12, 5.12, 5.12, 5.12 }, 1e-6, 100000, make_shared<f3_Mod3rdDeJong>() },
		//{ "f4_Mod4thDeJong_quartic", 10, 0.9, 0, f4_Mod4thDeJong_quartic::lower(), f4_Mod4thDeJong_quartic::upper(), 0, 100000, make_shared<f4_Mod4thDeJong_quartic>(rng) },
		{ "f5_Foxholes", 15, 0.9, 0, { -65.536, -65.536}, { 65.536, 65.536 }, 0.998005, 100000, make_shared<f5_Foxholes>() },
		{ "f6_Corana", 10, 0.5, 0, { -1000, -1000, -1000, -1000}, { 1000, 1000, 1000, 1000}, 1e-6, 100000, make_shared<f6_Corana>() },
		{ "f7_Griewangk", 25, 0.5, 0.2, vector<double>(10, -400), vector<double>(10, 400), 1e-6, 100000, make_shared<f7_Griewangk>() },
		{ "f8_Zimmermann", 10, 0.9, 0.9, { 0.0, 0.0 }, { 100.0, 100.0 }, 1e-6, 100000, make_shared<f8_Zimmermann>() },
		{ "f9_k4_Poly", 60, 0.6, 1, vector<double>(9,-100), vector<double>(9,100), 1e-6, 100000, make_shared<f9_k4_Poly>() },
		{ "f9_k8_Poly", 100, 0.6, 1, vector<double>(17,-1000), vector<double>(17,1000), 1e-6, 100000, make_shared<f9_k8_Poly>() },
		{ "f11_HyperEllipsoid_30", 20, 0.5, 0.1, vector<double>(30, -1), vector<double>(30, 1), 1e-10, 100000, make_shared<f11_HyperEllipsoid>(30) },
		{ "f11_HyperEllipsoid_100", 20, 0.5, 0.1, vector<double>(100, -1), vector<double>(100, 1), 1e-10, 100000, make_shared<f11_HyperEllipsoid>(100) },
		{ "f12_Katsuura_10", 15, 0.5, 0.1, vector<double>(10, -1000), vector<double>(10, 1000), 1.05, 100000, make_shared<f12_Katsuura>(10) },
		{ "f12_Katsuura_30", 15, 0.5, 0.1, vector<double>(30, -1000), vector<double>(30, 1000), 1.05, 100000, make_shared<f12_Katsuura>(30) },
	};

	auto comp = make_shared<ValueFitnessComparison<double>>();

	for (const auto &test : tests)
	{
		cout << test.name << ": "; 

		size_t successCount = 0;
		size_t totalCalc = 0;

		for (size_t run = 0 ; run < numRuns ; run++)
		{
			auto mut = make_shared<VectorDifferentialEvolutionMutation<double>>();
			auto cross = make_shared<VectorDifferentialEvolutionCrossover<double>>(rng);
			
			VectorDifferentialEvolutionIndividualCreation<double,double> creation(test.bottom, test.top, rng);
			ValueToReachStop<double> stop(test.VTR, test.maxGenerations);
		
			DifferentialEvolutionEvolver evolver(rng, mut, test.F, cross, test.CR, comp);
			SingleThreadedPopulationFitnessCalculation popCalc(test.calculator);

			MyEA ea;
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

				successCount++;
				totalCalc += test.calculator->getNumberOfEvaluations();
			}

			test.calculator->resetEvaluationCount();
		}

		cout << successCount << "/" << numRuns;
		cout << ", nfe " << totalCalc/numRuns << ", F = " << test.F << ", CR = " << test.CR << ", NP = " << test.popSize;
		if (successCount != numRuns)
			cout << " !!";
		cout << endl;
	}

	return 0;
}

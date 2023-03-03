#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cassert>
#include <limits>
#include <algorithm>
#include <tuple>
#include <memory>
#include <string>

using namespace std;

struct Problem
{
	virtual double evaluate(const vector<double> &trial) = 0;
};

template<size_t D>
struct f1_Sphere : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sumSquared = 0;
		for (auto v : x)
			sumSquared += v*v;
		return sumSquared;
	}
};

template<size_t D>
struct f2 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0, prod = 1;
		for (auto v : x)
		{
			sum += std::abs(v);
			prod *= std::abs(v);
		}
		return sum + prod;
	}
};

template<size_t D>
struct f3 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0;
		
		for (size_t i = 0 ; i < D ; i++)
		{
			double sub = 0;
			for (size_t j = 0 ; j <= i ; j++)
				sub += x[j];
			sum += sub*sub;
		}
		return sum;
	}
};

template<size_t D>
struct f4 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double max = -numeric_limits<double>::max();
		for (auto v : x)
			if (std::abs(v) > max)
				max = std::abs(v);
		return max;
	}
};

template<size_t D>
struct f5 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0;
		for (size_t i = 0 ; i < D-1 ; i++)
			sum += 100.0*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i]) + (x[i] - 1.0)*(x[i] - 1.0);
		return sum;
	}
};

template<size_t D>
struct f6 : public Problem
{
public:
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0;
		for (auto v : x)
			sum += std::floor(v+0.5)*std::floor(v+0.5);
		return sum;
	}
};

template<size_t D>
struct f7 : public Problem // Doesn't seem to work well
{
	f7(mt19937 &rng) : m_rng(rng) { }
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0;
		for (size_t i = 0 ; i < D ; i++)
			sum += (i+1.0)*x[i]*x[i]*x[i]*x[i];
		return sum + (uniform_real_distribution<>(0,1))(m_rng);
	}
private:
	mt19937 &m_rng;
};

template<size_t D>
struct f8 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);

		double sum = 0;
		for (auto v : x)
			sum += -v*std::sin(std::sqrt(std::abs(v)));
		return sum + 418.98288727243369*D;
	}
};

template<size_t D>
struct f9 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0;
		for (auto v : x)
			sum += v*v -10.0*std::cos(2.0*M_PI*v) + 10.0;
		return sum;
	}
};

template<size_t D>
struct f10 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum1 = 0, sum2 = 0;
		for (auto v : x)
		{
			sum1 += v*v;
			sum2 += std::cos(2.0*M_PI*v);
		}

		return -20.0*std::exp(-0.2*std::sqrt(sum1/D)) - std::exp(sum2/D) + 20 + std::exp(1);
	}
};

template<size_t D>
struct f11 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		double sum = 0, prod = 1;
		for (size_t i = 0 ; i < D ; i++)
		{
			sum += x[i]*x[i];
			prod *= std::cos(x[i]/std::sqrt(i+1.0));
		}
		return sum/4000.0 - prod + 1.0;
	}
};

inline double u(double z, double a, double k, double m)
{
	if (z > a)
		return k*std::pow(z-a, m);
	if (z < -a)
		return k*std::pow(-z-a, m);
	return 0;
};

template<size_t D>
struct f12 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);
		auto y = [&x](size_t i) { return 1.0+0.25*(x[i] + 1.0); };

		double sum = 10.0*std::sin(M_PI*y(0))*std::sin(M_PI*y(0));

		for (size_t i = 0 ; i < D-1 ; i++)
			sum += (y(i)-1.0)*(y(i)-1.0)*(1.0+10.0*std::sin(M_PI*y(i+1))*std::sin(M_PI*y(i+1)));

		sum += (y(D-1)-1.0)*(y(D-1)-1.0);
		sum *= M_PI/D;

		for (size_t i = 0 ; i < D ; i++)
			sum += u(x[i], 10, 100, 4);
		return sum;
	}
};

template<size_t D>
struct f13 : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == D);

		double sum = std::sin(3*M_PI*x[0])*std::sin(3*M_PI*x[0]);

		for (size_t i = 0 ; i < D-1 ; i++)
			sum += (x[i] - 1.0)*(x[i] - 1.0)*(1.0+std::sin(3*M_PI*x[i+1])*std::sin(3*M_PI*x[i+1]));

		sum += (x[D-1]-1.0)*(x[D-1]-1.0)*(1.0*std::sin(2*M_PI*x[D-1])*std::sin(2*M_PI*x[D-1]));
		sum *= 0.1;

		for (size_t i = 0 ; i < D ; i++)
			sum += u(x[i], 5, 100, 4);

		return sum;
	}
};

struct f14_Branin : public Problem
{
	double evaluate(const vector<double> &x) override
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

struct f15_GoldsteinPrice : public Problem
{
	double evaluate(const vector<double> &x) override
	{
		assert(x.size() == 2);

		return (1.0+std::pow(x[0]+x[1]+1.0,2)*(19.0-14.0*x[0]+3.0*x[0]*x[0]-14.0*x[1]+6.0*x[0]*x[1]+3.0*x[1]*x[1]))*(
			30.0+std::pow(2.0*x[0]-3.0*x[1],2)*(18.0-32.0*x[0]+12.0*x[0]*x[0]+48.0*x[1]-36.0*x[0]*x[1]+27.0*x[1]*x[1])
		);
	}
};

struct f16_Hartman3D : public Problem
{
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

	double evaluate(const vector<double> &x) override
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

struct f17_Hartman6D : public Problem
{
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

	double evaluate(const vector<double> &x) override
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

struct f_Shekel : public Problem
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

	double evaluate(const vector<double> &x) override
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

inline vector<double> generateIndividual(mt19937 &rng, const vector<pair<double,double>> &IPR)
{
	vector<double> x(IPR.size());
	for (size_t i = 0 ; i < x.size(); i++)
	{
		assert(IPR[i].first < IPR[i].second);
		x[i] = (uniform_real_distribution<>(IPR[i].first, IPR[i].second))(rng);
	}
	
	return x;
}

inline double clip(double x, double x0, double x1)
{
	if (x < x0)
		return x0;
	if (x > x1)
		return x1;
	return x;
}

inline double meanA(const vector<double> &P)
{
	if (P.size() == 0)
		return 0;

	double sum = 0;
	for (auto x : P)
		sum += x;
	sum /= P.size();
	return sum;
}

inline double meanL(const vector<double> &P)
{
	if (P.size() == 0)
		return 0;

	double sum = 0, sum2 = 0;
	for (auto x : P)
	{
		sum += x;
		sum2 += x*x;
	}
	return sum2/sum;
}

tuple<bool, vector<double>, double, size_t, size_t, vector<pair<double,double>>>
jade(mt19937 &rng, Problem &problem, size_t NP, const vector<pair<double,double>> &IPR,
         double VTR, size_t maxGen, bool useArchive, bool enforceBounds = false,
		 double p = 0.05, double c = 0.1)
{
	assert(IPR.size() > 0);
	assert(NP > 4);
	assert(p >= 0 && p <= 1);
	
	bool dumpPop = false;
	if (getenv("DUMPPOP"))
		dumpPop = true;

	size_t evaluations = 0;
	double muF = 0.5, muCR = 0.5;
	vector<pair<vector<double>,double>> P(NP), Q(NP);
	vector<pair<double,double>> muLog;

	for (auto &ind : P)
	{
		ind.first = generateIndividual(rng, IPR);
		ind.second = problem.evaluate(ind.first);
		evaluations++;
	}
	sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });

	auto logPop = [&P, dumpPop](size_t gen)
	{
		if (!dumpPop)
			return;

		cerr << "Generation " << gen << endl;
		for (size_t i = 0 ; i < P.size() ; i++)
		{
			cerr << "(" << i << ") [";
			for (auto x : P[i].first)
				cerr << " " << x;
			cerr << " ] => " << P[i].second << endl;
		}
	};

	logPop(0);

	vector<vector<double>> A;
	vector<double> SF, SCR;

	for (size_t g = 0 ; g < maxGen ; g++)
	{
		SCR.clear();
		SF.clear();

		muLog.push_back({ muF, muCR});
		for (size_t i = 0 ; i < NP ; i++)
		{
			normal_distribution<> randn(muCR, 0.1);
			cauchy_distribution<> randc(muF, 0.1);
			double CRi = clip(randn(rng), 0, 1);
			double Fi = clip(randc(rng), 0, 1);

			size_t xBest = (size_t)((uniform_real_distribution<>(0, p))(rng) * NP);
			size_t xR1 = 0, xR2 = 0;
			do { xR1 = (uniform_int_distribution<>(0, NP-1))(rng); } while (xR1 == i);
			do { xR2 = (uniform_int_distribution<>(0, NP+A.size()-1))(rng); } while (xR2 == i || xR2 == xR1);

			vector<double> v(IPR.size());
			auto &xi = P[i].first;
			if (xR2 < NP) // choose r2 from main population
			{
				for(size_t j = 0 ; j < v.size() ; j++)
				{
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - P[xR2].first[j]);
					assert(!isnan(v[j]));
				}
			}
			else // choose r2 from archive
			{
				assert(xR2 - NP < A.size());
				for(size_t j = 0 ; j < v.size() ; j++)
				{
					v[j] = xi[j] + Fi*(P[xBest].first[j] - xi[j]) + Fi*(P[xR1].first[j] - A[xR2-NP][j]);
					assert(!isnan(v[j]));
				}
			}

			vector<double> u(v.size());
			size_t jrand = (uniform_int_distribution<>(0, u.size()-1))(rng);
			for (size_t j = 0 ; j < u.size() ; j++)
			{
				if (j == jrand || (uniform_real_distribution<>(0,1))(rng) < CRi)
					u[j] = v[j];
				else
					u[j] = xi[j];
			}

			if (enforceBounds)
			{
				for (size_t j = 0 ; j < u.size() ; j++)
				{
					double xlow = IPR[j].first;
					double xhigh = IPR[j].second;
					if (u[j] < xlow)
						u[j] = (xlow + xi[j])/2.0;
					if (u[j] > xhigh)
						u[j] = (xhigh + xi[j])/2.0;
				}
			}

			double score = problem.evaluate(u);
			evaluations++;
			if (score < VTR)
			{
				return { true, u, score, g, evaluations, muLog };
			}

			if (P[i].second <= score)
			{
				// Nothing to do, just keep individual
				Q[i] = P[i];
			}
			else
			{
				if (useArchive)
					A.push_back(xi);

				Q[i].second = score;
				Q[i].first = u;
				SF.push_back(Fi);
				SCR.push_back(CRi);
			}
		}
		swap(P, Q);

		// Randomly remove solutions from A to trim size
		while(A.size() > NP)
		{
			size_t idx = (uniform_int_distribution<>(0, A.size()-1))(rng);
			if (idx != A.size()-1) // move last one to selected position
				A[idx] = A[A.size()-1];
			A.resize(A.size()-1);
		}

		// What to do if these sets are empty? Use 0 as mean? Don't update?
		muCR = (1.0-c)*muCR + c*meanA(SCR);
		muF = (1.0-c)*muF + c*meanL(SF);

		sort(P.begin(), P.end(), [](auto &i, auto &j) { return i.second < j.second; });
		logPop(g+1);
	}

	return { false, P[0].first, P[0].second, maxGen, evaluations, muLog };
}

struct Test
{
	string name;
	shared_ptr<Problem> problem;
	vector<pair<double,double>> IPR;
	double VTR;
	size_t NP;
	size_t maxGen;
	bool enforceBounds;
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

int main(int argc, char const *argv[])
{
	random_device rd;
	unsigned int seed = rd();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	mt19937 rng(seed);

	string matchName;
	if (argc > 1)
		matchName = argv[1];

	bool archive = true;
	if (getenv("NOARCH"))
		archive = false;

	vector<Test> tests {
		{ "f1_Sphere_30", make_shared<f1_Sphere<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f1_Sphere_100", make_shared<f1_Sphere<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
		{ "f2_30", make_shared<f2<30>>(), vector<pair<double,double>>(30,{-10,10}), 1e-8, 100, 100000, false },
		{ "f2_100", make_shared<f2<100>>(), vector<pair<double,double>>(100,{-10,10}), 1e-8, 400, 100000, false },
		{ "f3_30", make_shared<f3<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f3_100", make_shared<f3<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
		{ "f4_30", make_shared<f4<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f4_100", make_shared<f4<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
		{ "f5_30", make_shared<f5<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f5_100", make_shared<f5<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
		{ "f6_30", make_shared<f6<30>>(), vector<pair<double,double>>(30,{-100,100}), 1e-8, 100, 100000, false },
		{ "f6_100", make_shared<f6<100>>(), vector<pair<double,double>>(100,{-100,100}), 1e-8, 400, 100000, false },
		// { "f7_30", make_shared<f7<30>>(rng), vector<pair<double,double>>(30,{-1.28,1.28}), 1e-2, 100, 100000, false },
		// { "f7_100", make_shared<f7<100>>(rng), vector<pair<double,double>>(100,{-1.28,1.28}), 1e-2, 400, 100000, false },
		{ "f8_30", make_shared<f8<30>>(), vector<pair<double,double>>(30,{-500,500}), 1e-8, 100, 100000, true },
		{ "f8_100", make_shared<f8<100>>(), vector<pair<double,double>>(100,{-500,500}), 1e-8, 400, 100000, true },
		{ "f9_30", make_shared<f9<30>>(), vector<pair<double,double>>(30,{-5.12,5.12}), 1e-8, 100, 100000, false },
		{ "f9_100", make_shared<f9<100>>(), vector<pair<double,double>>(100,{-5.12,5.12}), 1e-8, 400, 100000, false },
		{ "f10_30", make_shared<f10<30>>(), vector<pair<double,double>>(30,{-32,32}), 1e-8, 100, 100000, false },
		{ "f10_100", make_shared<f10<100>>(), vector<pair<double,double>>(100,{-32,32}), 1e-8, 400, 100000, false },
		{ "f11_30", make_shared<f11<30>>(), vector<pair<double,double>>(30,{-600,600}), 1e-8, 100, 100000, false },
		{ "f11_100", make_shared<f11<100>>(), vector<pair<double,double>>(100,{-600,600}), 1e-8, 400, 100000, false },
		{ "f12_30", make_shared<f12<30>>(), vector<pair<double,double>>(30,{-50,50}), 1e-8, 100, 100000, false },
		{ "f12_100", make_shared<f12<100>>(), vector<pair<double,double>>(100,{-50,50}), 1e-8, 400, 100000, false },
		{ "f13_30", make_shared<f13<30>>(), vector<pair<double,double>>(30,{-50,50}), 1e-8, 100, 100000, false },
		{ "f13_100", make_shared<f13<100>>(), vector<pair<double,double>>(100,{-50,50}), 1e-8, 400, 100000, false },
		{ "f14_Branin", make_shared<f14_Branin>(), { {-5, 10}, {0, 15} }, 0.3978873577 + 1e-8, 30, 100000, false },
		{ "f15_GoldsteinPrice", make_shared<f15_GoldsteinPrice>(), { { -2, 2 }, { -2, 2 } }, 3.0 + 1e-8, 30, 100000, false },
		{ "f16_Hartman3D", make_shared<f16_Hartman3D>(), {{0,1},{0,1},{0,1}}, -3.862779787 + 1e-8, 30, 100000, false },
		{ "f17_Hartman6D", make_shared<f17_Hartman6D>(), {{0,1},{0,1},{0,1},{0,1},{0,1},{0,1}}, -3.322368011 + 1e-8, 30, 100000, false },
		{ "f18_Shekel5", make_shared<f_Shekel>(5), {{0,10},{0,10},{0,10},{0,10}}, -10.153199679 + 1e-8, 30, 100000, false },
		{ "f19_Shekel7", make_shared<f_Shekel>(7), {{0,10},{0,10},{0,10},{0,10}}, -10.402915336 + 1e-8, 30, 100000, false },
		{ "f20_Shekel10", make_shared<f_Shekel>(10), {{0,10},{0,10},{0,10},{0,10}}, -10.536443153 + 1e-8, 30, 100000, false },
	};
	
	size_t numRuns = 20;
	if (getenv("NUMRUNS"))
		numRuns = stoi(getenv("NUMRUNS"));

	bool printBest = false;
	if (getenv("DUMPFINAL"))
		printBest = true;

	size_t testCount = 0;
	for (auto &test : tests)
	{
		if (matchName.length() > 0 && matchName != test.name)
			continue;

		testCount++;

		AvgStd nfe;
		AvgStd genAvg;

		cout << test.name << ": ";
		for (size_t r = 0 ; r < numRuns ; r++)
		{
			bool converged;
			vector<double> best;
			double score;
			size_t generations, evaluations;
			vector<pair<double,double>> muLog;

			tie(converged, best, score, generations, evaluations, muLog) = jade(rng, *test.problem,
				test.NP, test.IPR, test.VTR, test.maxGen, archive, test.enforceBounds);

			if (converged)
			{
				nfe.process(evaluations);
				genAvg.process(generations);
				if (printBest)
				{
					cout << "FOUND: [";
					for (auto x : best)
						cout << " " << x;
					cout << " ] => " << score << endl;
				}
			}

			if (getenv("LOGFCR"))
				for (auto mu : muLog)
					cerr << mu.first << " " << mu.second << endl;
		}

		cout << " gen " << genAvg.toString();
		cout << " NFE " << nfe.toString();
		cout << " " << nfe.getCount() << "/" << numRuns;
		if (archive)
			cout << " archive";
		else
			cout << " no archive";

		cout << endl;
	}

	if (testCount == 0)
	{
		cerr << "No tests run, known are:" << endl;
		for (const auto &test : tests)
			cerr << "  " << test.name << endl;
	}
	return 0;
}

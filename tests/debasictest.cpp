#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <limits>
#include <cmath>
#include <map>
#include <memory>
#include <algorithm>
#include <string>
#include "testfunctions.h"

using namespace std;

struct Problem
{
	virtual size_t NP() = 0;
	virtual size_t D() = 0;
	virtual vector<vector<double>> IPR() = 0;
	virtual double F() = 0;
	virtual double CR() = 0;
	virtual double VTR() = 0;
	virtual double evaluate(const vector<double> &trial) = 0;
};

template <class T>
struct TestProblemTemplate : public Problem
{
	template <typename... Args>
	TestProblemTemplate(size_t NP, double F, double CR, double VTR, Args&&... args)
		: m_NP(NP), m_F(F), m_CR(CR), m_VTR(VTR), m_f(std::forward<Args>(args)...) { }
	vector<vector<double>> IPR() override
	{
		auto [ lower, upper ] = m_f.getInitialParameterRange();
		if (lower.size() != upper.size())
			throw runtime_error("Internal error: IPR ranges of different length");

		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < lower.size() ; i++)
			ipr.push_back({lower[i], upper[i]});

		return ipr;
	}

	size_t NP() override { return m_NP; }
	double F() override { return m_F; }
	double CR() override { return m_CR; }
	double VTR() override { return m_VTR; }

	size_t D() override { return m_f.getDimensions(); }
	double evaluate(const vector<double> &trial) override { return m_f.calculate1(trial); }

	T m_f;
	size_t m_NP;
	double m_F, m_CR, m_VTR;
};

typedef TestProblemTemplate<eatk::testfunctions::Sphere> f1_Sphere;
typedef TestProblemTemplate<eatk::testfunctions::Rosenbrock> f2_Rosenbrock;
typedef TestProblemTemplate<eatk::testfunctions::Mod3rdDeJong> f3_Step;
typedef TestProblemTemplate<eatk::testfunctions::QuarticWithNoise> f4_Quartic;
typedef TestProblemTemplate<eatk::testfunctions::Foxholes> f5_Foxholes;
typedef TestProblemTemplate<eatk::testfunctions::Corana> f6_Corana;
typedef TestProblemTemplate<eatk::testfunctions::Griewank> f7_Griewangk;
typedef TestProblemTemplate<eatk::testfunctions::Zimmermann> f8_Zimmermann;
typedef TestProblemTemplate<eatk::testfunctions::k4_Poly> f9_k4_Poly;
typedef TestProblemTemplate<eatk::testfunctions::k8_Poly> f9_k8_Poly;
typedef TestProblemTemplate<eatk::testfunctions::HyperEllipsoid> f11_HyperEllipsoid;
typedef TestProblemTemplate<eatk::testfunctions::Katsuura> f12_Katsuura;
typedef TestProblemTemplate<eatk::testfunctions::Rastrigin> f13_Rastrigin;
typedef TestProblemTemplate<eatk::testfunctions::Griewank> f14_Griewangk;
typedef TestProblemTemplate<eatk::testfunctions::AckleyFunction1> f15_Ackley;

struct f16_Goldstein : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 1; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 7.0 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		double X2 = X*X;
		double X4 = X2*X2;
		double X6 = X4*X2;
		return X6 - 15.0*X4 + 27.0*X2 + 250.0;
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

struct f17_PenalizedShubert : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 1; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return -12.8708855 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		
		double gX = 0;
		for (size_t i = 1 ; i <= 5 ; i++)
			gX += std::cos(X*(i+1)+i)*i;

		return gX+u(X, 10, 100, 2);
	}
};

struct f18_2DPenalizedShubert : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return -186.7309088 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double x1 = x[0];
		double x2 = x[1];
		
		double gX1 = 0;
		double gX2 = 0;
		for (size_t i = 1 ; i <= 5 ; i++)
		{
			gX1 += std::cos(x1*(i+1)+i)*i;
			gX2 += std::cos(x2*(i+1)+i)*i;
		}

		return gX1*gX2+u(x1, 10, 100, 2)+u(x2, 10, 100, 2);
	}
};

struct f19_Modified2DPenalizedShubert : public Problem
{
	const double m_beta;

	f19_Modified2DPenalizedShubert(double beta) : m_beta(beta) { }

	size_t NP() override { return 40; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 1.0; }
	double CR() override { return 0; }
	double VTR() override { return -186.7309088 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double x1 = x[0];
		double x2 = x[1];
		
		double gX1 = 0;
		double gX2 = 0;
		for (size_t i = 1 ; i <= 5 ; i++)
		{
			gX1 += std::cos(x1*(i+1)+i)*i;
			gX2 += std::cos(x2*(i+1)+i)*i;
		}

		double part1 = gX1*gX2+u(x1, 10, 100, 2)+u(x2, 10, 100, 2);
		return part1 + m_beta*((x1+1.42513)*(x1+1.42513) + (x2+0.80032)*(x2+0.80032));
	}
};

struct f20_SixHumpCamel : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return -1.0316285 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double x1 = x[0];
		double x2 = x[1];

		return (4.0 - 2.1*x1*x1 + x1*x1*x1*x1/3.0)*x1*x1 + x1*x2 + (-4.0+4.0*x2*x2)*x2*x2;
	}
};

struct f21_Function : public Problem
{
	const size_t m_D;

	f21_Function(size_t D) : m_D(D) { }
	
	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-10, 10});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 0 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double su = 0;
		for (double v : x)
			su += u(v, 10, 100, 4);

		double s = 10.0*std::pow(std::sin(M_PI+0.25*M_PI*(x[0]-1.0)),2);
		for (size_t i = 0 ; i < m_D-1 ; i++)
			s += 0.125*(x[i]-1.0)*(x[i]-1.0)*(1.0+10.0*std::pow(std::sin(M_PI+0.25*M_PI*(x[i+1]-1.0)),2));
		s += 0.125*(x[m_D-1]-1.0)*(x[m_D-1]-1.0);
		return s*M_PI/m_D + su;
	}
};

struct f22_Function : public Problem
{
	const size_t m_D;

	f22_Function(size_t D) : m_D(D) { }
	
	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-10, 10});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 0 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double su = 0;
		for (double v : x)
			su += u(v, 10, 100, 4);

		double s = 10.0*std::pow(std::sin(M_PI*x[0]),2);
		for (size_t i = 0 ; i < m_D-1 ; i++)
			s += (x[i]-1.0)*(x[i]-1.0)*(1.0+10.0*std::pow(std::sin(M_PI*x[i+1]),2));
		s += (x[m_D-1]-1.0)*(x[m_D-1]-1.0);
		return s*M_PI/m_D + su;
	}
};

struct f23_Function : public Problem
{
	const size_t m_D;

	f23_Function(size_t D) : m_D(D) { }
	
	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-10, 10});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 0 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double su = 0;
		for (double v : x)
			su += u(v, 10, 100, 4);

		double s = std::pow(std::sin(3.0*M_PI*x[0]),2);
		for (size_t i = 0 ; i < m_D-1 ; i++)
			s += (x[i]-1.0)*(x[i]-1.0)*(1.0+std::pow(std::sin(3.0*M_PI*x[i+1]),2));
		s += (x[m_D-1]-1.0)*(x[m_D-1]-1.0)*(1.0*std::pow(std::sin(2.0*M_PI*x[m_D-1]),2));
		
		return s*0.1 + su;
	}
};

struct f24_Function : public Problem
{
	const size_t m_D;

	f24_Function(size_t D) : m_D(D) { }
	
	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-10, 10});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 0 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double su = 0;
		for (double v : x)
			su += u(v, 5, 100, 4);

		double s = std::pow(std::sin(3.0*M_PI*x[0]),2);
		for (size_t i = 0 ; i < m_D-1 ; i++)
			s += (x[i]-1.0)*(x[i]-1.0)*(1.0+std::pow(std::sin(3.0*M_PI*x[i+1]),2));
		s += (x[m_D-1]-1.0)*(x[m_D-1]-1.0)*(1.0*std::pow(std::sin(2.0*M_PI*x[m_D-1]),2));
		
		return s*0.1 + su;
	}
};

struct f25_Function : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 1; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return -0.3523861 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		return 0.25*X*X*X*X -0.5*X*X + 0.1*X;
	}
};

struct f26_Function : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return -0.3523861 + 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		double Y = x[1];
		return 0.25*X*X*X*X -0.5*X*X + 0.1*X + 0.5*Y*Y;
	}
};

struct f27_Function : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		double Y = x[1];
		return 0.5*X*X + 0.5*(1.0-std::cos(2.0*X)) + Y*Y;
	}
};

struct f28_Function : public Problem
{
	const int m_m, m_n;

	f28_Function(int n, int m) : m_n(n), m_m(m) { }
	size_t NP() override { return 20; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override
	{ 
		double base = 0;
		if (m_n == 1 && m_m == -1)
			base = -0.4074616;
		else if (m_n == 2 && m_m == -2)
			base = -18.0586967;
		else if (m_n == 3 && m_m == -3)
			base = -227.7657500;
		else if (m_n == 4 && m_m == -4)
			base = -2429.4147670;
		else if (m_n == 5 && m_m == -5)
			base = -24776.5183423;
		else if (m_n == 6 && m_m == -6)
			base = -249293.0182630;
		else
			throw runtime_error("Invalid m and n settings");
		return base + 1e-6; 
	}
	double evaluate(const vector<double> &x) override
	{
		double X = x[0];
		double Y = x[1];
		return std::pow(10.0, m_n)*X*X + Y*Y - std::pow(X*X+Y*Y,2) + std::pow(10.0, m_m)*std::pow(X*X+Y*Y, 4);
	}
};

struct f29_Function : public Problem
{
	size_t NP() override { return 20; }
	size_t D() override { return 5; }
	vector<vector<double>> IPR() override
	{
		return { { -10, 10 }, { -10, 10 }, { -10, 10 }, { -10, 10 }, { -10, 10 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 1e-6; } 
	double evaluate(const vector<double> &x) override
	{
		double r = 0;
		for (size_t i = 0 ; i < 5 ; i++)
			r += (1.0+i)*x[i]*x[i];
		return std::pow(r, 0.25);
	}
};

struct f30_Function : public Problem // Doesn't seem to work?
{
	static const vector<double> z;
	static const vector<double> delta;

	size_t NP() override { return 30; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -std::pow(std::exp(1), 4), { std::pow(std::exp(1), 4) } }, { -std::pow(std::exp(1), 4), { std::pow(std::exp(1), 4) } } };
	}
	double F() override { return 0.5; }
	double CR() override { return 1.0; }
	double VTR() override { return -0.000888085 + 1e-6; } 
	double evaluate(const vector<double> &x) override
	{
		double Fx = 1.0;
		auto phi = [](double z)
		{
			return 0.5*(1.0 + std::erf(z/std::sqrt(2)));
		};
		for (size_t i = 0 ; i < 14 ; i++)
			Fx *= std::pow(phi(z[i]-x[0])/x[1],1.0-delta[i])*std::pow(1.0-phi(z[i]-x[0])/x[1],delta[i]);
		
		return -Fx + u(x[0], 10000, 100, 2) + u(x[1], 10000, 100, 2);
	}
};

const vector<double> f30_Function::z = { 1219, 1371, 1377, 1144, 1201, 1225, 1244, 1254, 1304, 1328, 1351, 1356, 1370, 1390 };
const vector<double> f30_Function::delta = { 0,0,0, 1,1,1, 1,1,1, 1,1,1, 1,1 };

bool runDEonProblem(mt19937 &rng, Problem &problem)
{
	uniform_real_distribution<> dist(0, 1);
	auto rnd_uni = [&rng, &dist]()
	{
		return dist(rng);
	};

	size_t count = 0, gen_max = 100000;
	size_t evaluation_count = 0;

	size_t NP = problem.NP();
	size_t D = problem.D();
	vector<vector<double>> IPR = problem.IPR();
	double F = problem.F();
	double CR = problem.CR();
	double VTR = problem.VTR();

	vector<vector<double>> x1(NP), x2(NP);
	vector<double> trial(D);
	vector<double> cost(NP);

	// Initialize the population vectors
	for (size_t i = 0 ; i < NP ; i++)
	{
		for (size_t j = 0 ; j < D ; j++)
		{
			x1[i].resize(D);
			x2[i].resize(D);
			x1[i][j] = (rnd_uni()*(IPR[j][1] - IPR[j][0]) + IPR[j][0]);
		}

		cost[i] = numeric_limits<double>::max(); 
	}

	while (count < gen_max)
	{
		for (size_t i = 0 ; i < NP ; i++)
		{
			size_t a, b, c;
			do { a = (size_t)(rnd_uni()*NP); } while (a == i);
			do { b = (size_t)(rnd_uni()*NP); } while (b == i || b == a );
			do { c = (size_t)(rnd_uni()*NP); } while (c == i || c == a || c == b);
			
			size_t j = (size_t)(rnd_uni()*D);
			for (size_t k = 1 ; k <= D ; k++)
			{
				if (rnd_uni() < CR || k == D)
					trial[j] = x1[c][j] + F*(x1[a][j] - x1[b][j]);
				else
					trial[j] = x1[i][j];

				j = (j + 1)%D;
			}

			double score = problem.evaluate(trial);
			evaluation_count++;

			if (score < VTR)
			{
				cout << "NP: " << NP << " D: " << D << " F: " << F << " CR: " << CR;
				cout << " gen: " << count << " eval: " << evaluation_count;
				cout << " score: " << score;
				cout << " sol: [";
				for (auto v : trial)
					cout << " " << v;
				cout << " ]" << endl;
				return true;
			}

			if (score <= cost[i])
			{
				// for (auto x : trial)
				//	 cerr << " " << x;
				// cerr << " => " << score << endl;

				for (size_t j = 0 ; j < D ; j++)
					x2[i][j] = trial[j];
				cost[i] = score;
			}
			else
			{
				for (size_t j = 0 ; j < D ; j++)
					x2[i][j] = x1[i][j];
			}
		}

		swap(x1, x2);

		count++;
	}

	cout << "No solution found after " << gen_max << " generations" << endl;
	return -1;

}

int main(int argc, char *argv[])
{
	random_device rd;
	unsigned int seed = rd();
	if (getenv("SEED"))
		seed = (unsigned int)stoul(getenv("SEED"));

	auto rng = make_shared<mt19937>(seed);
	
	map<string, shared_ptr<Problem>> problems {
		{ "f1", make_shared<f1_Sphere>(5, 0.9, 0.1, 1e-6, 3, pair(-5.12, 5.12)) },
		{ "f2", make_shared<f2_Rosenbrock>(10, 0.9, 0.9, 1e-6, pair(-2.048, 2.048)) },
		{ "f3", make_shared<f3_Step>(10, 0.9, 0.0, 1e-6, pair(-5.12, 5.12)) },
		//{ "f4", make_shared<f4_Quartic>(10, 0.9, 0.0, 15.0, 30, rng, true, pair(-1.28, 1.28)) },
		{ "f5", make_shared<f5_Foxholes>(15, 0.9, 0.0, 0.998005, pair(-65.536, 65.536)) },
		{ "f6", make_shared<f6_Corana>(10, 0.5, 0.0, 1e-6, pair(-1000.0, 1000.0)) },
		{ "f7", make_shared<f7_Griewangk>(25, 0.5, 0.2, 1e-6, 10, pair(-400.0,400.0)) },
		{ "f8", make_shared<f8_Zimmermann>(10, 0.9, 0.9, 1e-6, pair(0.0, 100.0)) },
		{ "f9_k4", make_shared<f9_k4_Poly>(60, 0.6, 1.0, 1e-6, pair(-100.0, 100.0)) },
		{ "f9_k8", make_shared<f9_k8_Poly>(100, 0.6, 1.0, 1e-6, pair(-1000.0, 1000.0)) },
		{ "f11_30", make_shared<f11_HyperEllipsoid>(20, 0.5, 0.1, 1e-10, 30, pair(-1.0, 1.0)) },
		{ "f11_100", make_shared<f11_HyperEllipsoid>(20, 0.5, 0.1, 1e-10, 100, pair(-1.0, 1.0)) },
		{ "f12_10", make_shared<f12_Katsuura>(15, 0.5, 0.1, 1.05, 10, pair(-1000.0, 1000.0)) },
		{ "f12_30", make_shared<f12_Katsuura>(15, 0.5, 0.1, 1.05, 30, pair(-1000.0, 1000.0)) },
		{ "f13_20", make_shared<f13_Rastrigin>(25, 0.5, 0.0, 0.9, 20, pair(-600.0, 600.0)) },
		{ "f13_100", make_shared<f13_Rastrigin>(25, 0.5, 0.0, 0.9, 100, pair(-600.0, 600.0)) },
		{ "f14_20", make_shared<f14_Griewangk>(20, 0.5, 0.1, 1e-3, 20, pair(-600.0, 600.0)) },
		{ "f14_100", make_shared<f14_Griewangk>(20, 0.5, 0.1, 1e-3, 100, pair(-600.0, 600.0)) },
		{ "f15_30", make_shared<f15_Ackley>(20, 0.5, 0.1, 1e-3, 30, -0.02, pair(-30.0, 30.0)) },
		{ "f15_100", make_shared<f15_Ackley>(20, 0.5, 0.1, 1e-3, 100, -0.02, pair(-30.0, 30.0)) },
		{ "f16", make_shared<f16_Goldstein>() },
		{ "f17", make_shared<f17_PenalizedShubert>() },
		{ "f18", make_shared<f18_2DPenalizedShubert>() },
		{ "f19_0.5", make_shared<f19_Modified2DPenalizedShubert>(0.5) },
		{ "f19_1", make_shared<f19_Modified2DPenalizedShubert>(1.0) },
		{ "f20", make_shared<f20_SixHumpCamel>() },
		{ "f21_2", make_shared<f21_Function>(2) },
		{ "f21_3", make_shared<f21_Function>(3) },
		{ "f21_4", make_shared<f21_Function>(4) },
		{ "f22_5", make_shared<f22_Function>(5) },
		{ "f22_8", make_shared<f22_Function>(8) },
		{ "f22_10", make_shared<f22_Function>(10) },
		{ "f23_2", make_shared<f23_Function>(2) },
		{ "f23_3", make_shared<f23_Function>(3) },
		{ "f23_4", make_shared<f23_Function>(4) },
		{ "f24_5", make_shared<f24_Function>(5) },
		{ "f24_6", make_shared<f24_Function>(6) },
		{ "f24_7", make_shared<f24_Function>(7) },
		{ "f25", make_shared<f25_Function>() },
		{ "f26", make_shared<f26_Function>() },
		{ "f27", make_shared<f27_Function>() },
		{ "f28_1", make_shared<f28_Function>(1, -1) },
		{ "f28_2", make_shared<f28_Function>(2, -2) },
		{ "f28_3", make_shared<f28_Function>(3, -3) },
		{ "f28_4", make_shared<f28_Function>(4, -4) },
		{ "f28_5", make_shared<f28_Function>(5, -5) },
		{ "f28_6", make_shared<f28_Function>(6, -6) },
		{ "f29", make_shared<f29_Function>() },
		//{ "f30", make_shared<f30_Function>() },
	};

	//f1_Sphere problem;
	//f2_Rosenbrock problem;
	//f3_Step problem;
	//f5_Foxholes problem;
	//f6_Corana problem;
	//f7_Griewangk problem;
	//f8_Zimmermann problem;
	//f9_k4_Poly problem;
	//f9_k8_Poly problem;
	//f11_HyperEllipsoid problem(30);
	//f11_HyperEllipsoid problem(100);
	//f12_Katsuura problem(10);
	//f12_Katsuura problem(30);
	//f13_Rastrigin problem(20);
	//f13_Rastrigin problem(100);
	//f14_Griewangk problem(20);
	//f14_Griewangk problem(100);
	//f15_Ackley problem(30);
	//f15_Ackley problem(100);
	//f16_Goldstein problem;
	//f17_PenalizedShubert problem;
	//f18_2DPenalizedShubert problem;
	//f19_Modified2DPenalizedShubert problem(0.5);
	//f19_Modified2DPenalizedShubert problem(1.0);
	//f20_SixHumpCamel problem;
	//f21_Function problem(2);
	//f21_Function problem(3);
	//f21_Function problem(4);
	//f22_Function problem(5);
	//f22_Function problem(8);
	//f22_Function problem(10);
	//f23_Function problem(2);
	//f23_Function problem(3);
	//f23_Function problem(4);
	//f24_Function problem(5);
	//f24_Function problem(6);
	//f24_Function problem(7);
	//f25_Function problem;
	//f26_Function problem;
	//f27_Function problem;
	//f28_Function problem(1, -1);
	//f28_Function problem(2, -2);
	//f28_Function problem(3, -3);
	//f28_Function problem(4, -4);
	//f28_Function problem(5, -5);
	//f28_Function problem(6, -6);
	//f29_Function problem;
	//f30_Function problem; // Doesn't seem to work well

	auto listProblems = [&problems]()
	{
		vector<string> names;
		for (auto it : problems)
			names.push_back(it.first);
		for (auto &n : names)
			cerr << "    " << n << endl;
	};

	if (argc != 2)
	{
		cerr << "Specify a problem name!" << endl;
		listProblems();
		return -1;
	}

	auto it = problems.find(argv[1]);
	if (it == problems.end())
	{
		cerr << "Specified problem name not found" << endl;
		listProblems();
		return -1;
	}

	if (!(runDEonProblem(*rng, *(it->second))))
		return -1;
	return 0;
}

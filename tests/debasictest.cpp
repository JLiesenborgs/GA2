#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <limits>
#include <cmath>
#include <map>
#include <memory>
#include <algorithm>
#include <string>

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

struct f1_Sphere : public Problem // Parameters from deshort1.ps
{
	size_t NP() override { return 5; }
	size_t D() override { return 3; }
	vector<vector<double>> IPR() override
	{
		return { { -5.12, 5.12 }, { -5.12, 5.12 }, { -5.12, 5.12 } };
	}
	double F() override { return 0.9; }
	double CR() override { return 0.1; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	}
};

struct f2_Rosenbrock : public Problem
{
	size_t NP() override { return 10; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -2.048, 2.048 }, { -2.048, 2.048 } };
	}
	double F() override { return 0.9; }
	double CR() override { return 0.9; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		return 100.0*(x[0]*x[0] - x[1])*(x[0]*x[0] - x[1]) + (1.0-x[0])*(1.0-x[0]);
	}
};

struct f3_Step : public Problem
{
	size_t NP() override { return 10; }
	size_t D() override { return 5; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < 5 ; i++)
			ipr.push_back({-5.12, 5.12});
		return ipr;
	}
	double F() override { return 0.9; }
	double CR() override { return 0.0; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double r = 30;
		for (double c: x)
		{
			if (std::abs(c) <= 5.12)
				r += std::floor(c);
			else if (c > 5.12)
				r += 30.0*(c-5.12);
			else
				r += 30.0*(5.12-c);
		}
		return r;
	}
};

struct f4_Quartic : public Problem // Is this right?
{
	f4_Quartic(mt19937 &rng) : m_rng(rng) { }

	size_t NP() override { return 10; }
	size_t D() override { return 30; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < 30 ; i++)
			ipr.push_back({-1.28, 1.28});
		return ipr;
	}
	double F() override { return 0.9; }
	double CR() override { return 0.0; }
	double VTR() override { return 15.0; }
	double evaluate(const vector<double> &x) override
	{
		uniform_real_distribution<> dist(0, 1);

		double s = 0;
		for (size_t j = 0 ; j < 30 ; j++)
		{
			double eta = dist(m_rng);
			s += (x[j]*x[j]*x[j]*x[j]*(j+1) + eta);
		}
		return s;

	}

	mt19937 &m_rng;
};

struct f5_Foxholes : public Problem
{
	vector<vector<double>> a;
	f5_Foxholes()
	{
		a = { 
			{ -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32,-32, -16, 0, 16, 32,-32, -16, 0, 16, 32},
			{ -32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0,0,0,0,0, 16,16,16,16,16, 32,32,32,32,32 }	  
		};
	}
	size_t NP() override { return 15; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { -65.536, 65.536 }, { -65.536, 65.536 } };
	}
	double F() override { return 0.9; }
	double CR() override { return 0; }
	double VTR() override { return 0.998005; }
	double evaluate(const vector<double> &x) override
	{
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
};

struct f6_Corana : public Problem
{
	size_t NP() override { return 10; }
	size_t D() override { return 4; }
	vector<vector<double>> IPR() override
	{
		return { { -1000, 1000 }, { -1000, 1000 }, { -1000, 1000 }, { -1000, 1000 } };
	}
	double F() override { return 0.5; }
	double CR() override { return 0; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
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

struct f7_Griewangk : public Problem
{
	size_t NP() override { return 25; }
	size_t D() override { return 10; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < 10 ; i++)
			ipr.push_back({-400,400});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0.2; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
		double s = 1.0;
		for (size_t i = 0 ; i < 10 ; i++)
			s+= x[i]*x[i]/4000.0;

		double p = 1.0;
		for (size_t i = 0 ; i < 10 ; i++)
			p *= std::cos(x[i]/std::sqrt(i+1));
		
		return s - p;
	}
};

struct f8_Zimmermann : public Problem // (settings from deshort1.ps) back to settings from paper
{
	size_t NP() override { return 10; }
	size_t D() override { return 2; }
	vector<vector<double>> IPR() override
	{
		return { { 0, 100 }, { 0, 100 } };
	}
	double F() override { return 0.9; }
	double CR() override { return 0.9; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
	{
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

struct f9_k4_Poly : public Problem
{
	vector<double> z;

	f9_k4_Poly()
	{
		z.push_back(-1.2);
		for (size_t i = 0 ; i < 60 ; i++)
			z.push_back(i*2.0/(60-1) + (-1.0));
		z.push_back(1.2);
	}
	size_t NP() override { return 60; }
	size_t D() override { return 9; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < 9 ; i ++)
			ipr.push_back({ -100, 100});
		return ipr;
	}
	double F() override { return 0.6; }
	double CR() override { return 1; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
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

struct f9_k8_Poly : public Problem
{
	vector<double> z;

	f9_k8_Poly()
	{
		z.push_back(-1.2);
		for (size_t i = 0 ; i < 60 ; i++)
			z.push_back(i*2.0/(60-1) + (-1.0));
		z.push_back(1.2);
	}
	size_t NP() override { return 100; }
	size_t D() override { return 17; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < 17 ; i ++)
			ipr.push_back({ -1000, 1000});
		return ipr;
	}
	double F() override { return 0.6; }
	double CR() override { return 1; }
	double VTR() override { return 1e-6; }
	double evaluate(const vector<double> &x) override
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

struct f11_HyperEllipsoid : public Problem
{
	const size_t m_D;

	f11_HyperEllipsoid(size_t D) : m_D(D) { }

	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i ++)
			ipr.push_back({ -1, 1 });
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0.1; }
	double VTR() override { return 1e-10; }
	double evaluate(const vector<double> &x) override
	{
		double s = 0;

		for (size_t j = 0 ; j < m_D ; j++)
			s += (j+1.0)*(j+1.0)*x[j]*x[j];

		return s;
	}
};

struct f12_Katsuura : public Problem
{
	const size_t m_D;

	f12_Katsuura(size_t D) : m_D(D) { }

	size_t NP() override { return 15; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i ++)
			ipr.push_back({ -1000, 1000 });
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0.1; }
	double VTR() override { return 1.05; }
	double evaluate(const vector<double> &x) override
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

struct f13_Rastrigin : public Problem
{
	const size_t m_D;

	f13_Rastrigin(size_t D) : m_D(D) { }

	size_t NP() override { return 25; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i ++)
			ipr.push_back({ -600, 600 });
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0.0; }
	double VTR() override { return 0.9; }
	double evaluate(const vector<double> &x) override
	{
		double s = 10.0*m_D;
		for (size_t j = 0 ; j < m_D ; j++)
			s += ( x[j]*x[j] - 10.0*std::cos(2.0*M_PI*x[j]) );

		return s;
	}
};

struct f14_Griewangk : public Problem
{
	const size_t m_D;

	f14_Griewangk(size_t D) : m_D(D) { }

	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-600,600});
		return ipr;
	}
	double F() override { return 0.5; }
	double CR() override { return 0.1; }
	double VTR() override { return 1e-3; }
	double evaluate(const vector<double> &x) override
	{
		double s = 1.0;
		for (size_t i = 0 ; i < m_D ; i++)
			s+= x[i]*x[i]/4000.0;

		double p = 1.0;
		for (size_t i = 0 ; i < m_D ; i++)
			p *= std::cos(x[i]/std::sqrt(i+1));
		
		return s - p;
	}
};

struct f15_Ackley : public Problem
{
	const size_t m_D;

	f15_Ackley(size_t D) : m_D(D) { }

	size_t NP() override { return 20; }
	size_t D() override { return m_D; }
	vector<vector<double>> IPR() override
	{
		vector<vector<double>> ipr;
		for (size_t i = 0 ; i < m_D ; i++)
			ipr.push_back({-30,30});
		return ipr;
	}
	double F() override { return 0.5; } // settings from article don't seem to work well?
	double CR() override { return 0.1; }
	double VTR() override { return 1e-3; }
	double evaluate(const vector<double> &x) override
	{
		double squaredSum = 0.0;
		double cosSum = 0.0;
		for (size_t i = 0 ; i < m_D ; i++)
		{
			squaredSum += x[i]*x[i];
			cosSum += std::cos(2.0*M_PI*x[i]);
		}

		return -20.0*std::exp(-0.02*std::sqrt(squaredSum/m_D)) - std::exp(cosSum/m_D) + 20.0 + std::exp(1.0);
	}
};

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
	mt19937 rng(rd());
	
	map<string, shared_ptr<Problem>> problems {
		{ "f1", make_shared<f1_Sphere>() },
		{ "f2", make_shared<f2_Rosenbrock>() },
		{ "f3", make_shared<f3_Step>() },
		// { "f4", make_shared<f4_Quartic>(rng) },
		{ "f5", make_shared<f5_Foxholes>() },
		{ "f6", make_shared<f6_Corana>() },
		{ "f7", make_shared<f7_Griewangk>() },
		{ "f8", make_shared<f8_Zimmermann>() },
		{ "f9_k4", make_shared<f9_k4_Poly>() },
		{ "f9_k8", make_shared<f9_k8_Poly>() },
		{ "f11_30", make_shared<f11_HyperEllipsoid>(30) },
		{ "f11_100", make_shared<f11_HyperEllipsoid>(100) },
		{ "f12_10", make_shared<f12_Katsuura>(10) },
		{ "f12_30", make_shared<f12_Katsuura>(30) },
		{ "f13_20", make_shared<f13_Rastrigin>(20) },
		{ "f13_100", make_shared<f13_Rastrigin>(100) },
		{ "f14_20", make_shared<f14_Griewangk>(20) },
		{ "f14_100", make_shared<f14_Griewangk>(100) },
		{ "f15_30", make_shared<f15_Ackley>(30) },
		{ "f15_100", make_shared<f15_Ackley>(100) },
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
		sort(names.begin(), names.end());
		cerr << "Available function names:" << endl;
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

	if (!(runDEonProblem(rng, *(it->second))))
		return -1;
	return 0;
}

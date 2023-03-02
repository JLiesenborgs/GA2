#include "mersennerandomnumbergenerator.h"
#include <random>
#include <iostream>

using namespace std;
using namespace eatk;

template<class Distribution>
inline double chooseTruncatedDistribution(RandomNumberGenerator &rng, double mu, double sigma, double x0, double x1)
{
	struct G
	{
		typedef uint64_t result_type;

		G(RandomNumberGenerator &_rng) : rng(_rng) { }
		uint64_t operator()() { return (uint64_t)rng.getRandomUint32(); }
		constexpr static uint64_t min() { return 0; }
		constexpr static uint64_t max() { return 0x100000000; }

		RandomNumberGenerator &rng;
	};

	G gen(rng);
	Distribution dist(mu, sigma);

	// int count = 0;
	// int maxCount = 100;
	double x = dist(gen);
	// while ((x > x1 || x < x0) && count < maxCount)
	// {
		// x = dist(gen);
		// count++;
	// }
	if (x < x0)
		x = x0;
	if (x > x1)
		x = x1;
	return x;
}

int main(int argc, char *argv[])
{
	random_device rd;
	MersenneRandomNumberGenerator rng(rd());
	//mt19937 rng(rd());

	if (argc != 2)
	{
		cerr << "Specify number of samples" << endl;
		return -1;
	}
	int samples = stoi(argv[1]);

	for (int i = 0 ; i < samples ; i++)
	{
		double x = chooseTruncatedDistribution<cauchy_distribution<>>(rng, 5, 2, -1, 10);
		cerr << x << endl;
	}

	return 0;
}

#include <cassert>
#include <iostream>
#include <string>
#include "mersennerandomnumbergenerator.h"

using namespace std;
using namespace eatk;

int main(int argc, char const *argv[])
{
    MersenneRandomNumberGenerator m_rng(12345);
	auto pickIndex = [&m_rng](size_t num) { return ((size_t)m_rng.getRandomUint32())%(num); };
	auto adjust = [](size_t &idx, size_t refIdx) { if (idx >= refIdx) idx++; };

	size_t targetPopulationSize = 16;

	auto pickRandomIndices = [targetPopulationSize,pickIndex,adjust](size_t i, size_t &r1, size_t &r2, size_t &r3)
	{
		r1 = pickIndex(targetPopulationSize-1);
		r2 = pickIndex(targetPopulationSize-2);
		r3 = pickIndex(targetPopulationSize-3);
		adjust(r3, r2);

		adjust(r2, r1);
		adjust(r3, r1);

		adjust(r1, i);
		adjust(r2, i);
		adjust(r3, i);

		assert(r1 != i && r2 != i && r3 != i);
		assert(r2 != r1 && r3 != r1);
		assert(r3 != r2);
	};

	size_t i = stoi(argv[1]);
	vector<size_t> countR1(targetPopulationSize, 0), countR2(targetPopulationSize, 0), countR3(targetPopulationSize, 0);
	for (size_t x = 0 ; x < 10000000 ; x++)
	{
		size_t r1, r2, r3;
		pickRandomIndices(i, r1, r2, r3);
		countR1[r1]++;
		countR2[r2]++;
		countR3[r3]++;
	}

	cout << "# i = 5" << endl;
	
	auto dump = [](const vector<size_t> &bins)
	{
		for (size_t i = 0 ; i < bins.size() ; i++)
		{
			cout << i << " " << bins[i] << endl;
			cout << (i+1) << " " << bins[i] << endl;
		}
		cout << endl;
		cout << endl;
	};

	dump(countR1);
	dump(countR2);
	dump(countR3);

	return 0;
}

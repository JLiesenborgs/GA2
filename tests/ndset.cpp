#include "population.h"
#include "mersennerandomnumbergenerator.h"
#include "vectorgenomefitness.h"
#include "basicnondominatedsetcreator.h"
#include "fasternondominatedsetcreator.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <assert.h>

using namespace std;
using namespace errut;
using namespace eatk;

int main(int argc, char const *argv[])
{
	auto error = [](const string &s)
	{
		cerr << s << endl;
		exit(-1);
	};

	random_device rd;
	MersenneRandomNumberGenerator rng(rd());

	Population pop, pop2;
	size_t N = 50;
	size_t dim = 2;

	auto buildPopulation = [N, dim, &rng](Population &pop)
	{
		for (size_t i = 0 ; i < N ; i++)
		{
			VectorGenome<float> dummy;
			VectorFitness<float> f { dim };
			f.setCalculated();

			for (size_t j = 0 ; j < dim ; j++)
				f.getValues()[j] = (int)(rng.getRandomDouble()*100);

			pop.append(make_shared<Individual>(dummy.createCopy(), f.createCopy()));
		}
	};
	buildPopulation(pop);
	buildPopulation(pop2); // To check merging of ND sets

	vector<shared_ptr<Individual>> ndSet, ndSet2;
	vector<shared_ptr<Individual>> remainder, remainder2;

	BasicNonDominatedSetCreator bndsc(make_shared<VectorFitnessComparison<float>>(), dim);
	FasterNonDominatedSetCreator fndsc(make_shared<VectorFitnessComparison<float>>(), dim);
	bool_t r;

	auto sortSet = [](auto &s)
	{
		sort(s.begin(), s.end(), [](auto &e1, auto &e2) {
			const VectorFitness<float> &f1 = dynamic_cast<const VectorFitness<float> &>(e1->fitnessRef());
			const VectorFitness<float> &f2 = dynamic_cast<const VectorFitness<float> &>(e2->fitnessRef());
			return f1.getValues()[0] < f2.getValues()[0];
		});
	};

	auto printSet = [dim](auto &v)
	{
		for (auto &i : v)
		{
			const VectorFitness<float> &f = dynamic_cast<const VectorFitness<float> &>(i->fitnessRef());
			for (size_t k = 0 ; k < dim ; k++)
				cout << f.getValues()[k] << " ";
			cout << endl;
		};
	};

	if (!(r = bndsc.calculateNonDomitatedSet(pop.individuals(), ndSet, remainder)))
		 error("Error: " + r.getErrorString());
	if (!(r = fndsc.calculateNonDomitatedSet(pop.individuals(), ndSet2, remainder2)))
		 error("Error: " + r.getErrorString());

	cout << "# ND set size: " << ndSet.size() << ", " << ndSet2.size() << endl;
	assert(ndSet.size() == ndSet2.size());
	assert(remainder.size() == remainder2.size());
	
	sortSet(ndSet);
	sortSet(ndSet2);
	sortSet(remainder);
	sortSet(remainder2);
	
	for (size_t i = 0 ; i < ndSet.size() ; i++)
		assert(ndSet[i].get() == ndSet2[i].get());

	for (size_t i = 0 ; i < remainder.size() ; i++)
		assert(remainder[i].get() == remainder2[i].get());

	printSet(pop.individuals());
	cout << endl << endl;
	// printSet(ndSet);
	// cout << endl << endl;
	// printSet(remainder);
	// cout << endl << endl;

	if (!(r = bndsc.calculateAllNDSets(pop.individuals())))
		error("Error: " + r.getErrorString());

	if (!(r = fndsc.calculateAllNDSets(pop.individuals())))
		error("Error: " + r.getErrorString());

	cout << "# Found " << bndsc.getNumberOfSets() <<  " sets" << endl;
	cout << "# Found " << fndsc.getNumberOfSets() <<  " sets" << endl;

	if (bndsc.getNumberOfSets() != fndsc.getNumberOfSets())
		error("Number of ND sets is not equal");

	for (size_t i = 0 ; i < bndsc.getNumberOfSets() ; i++)
	{
		auto s = bndsc.getSet(i);
		auto s2 = fndsc.getSet(i);

		sortSet(s);
		sortSet(s2);

		assert(s.size() == s2.size());
		for (size_t j = 0 ; j < s.size() ; j++)
			assert(s[j].get() == s2[j].get());

		printSet(s);
		cout << endl << endl;
	}

	if (!(r = bndsc.calculateNonDomitatedSet(pop.individuals(), ndSet2, remainder2)))
		 error("Error: " + r.getErrorString());

	vector<shared_ptr<Individual>> merged = ndSet, merged2 = ndSet;

	if (!(r = bndsc.mergeNDSets(merged, ndSet2)))
		error("Can't merge: " + r.getErrorString());
	
	if (!(r = fndsc.mergeNDSets(merged2, ndSet2)))
		error("Can't merge: " + r.getErrorString());

	cout << "# merged ND sets sizes: " << merged.size() << ", " << merged2.size() << endl;
	assert(merged.size() == merged2.size());
	for (size_t i = 0 ; i < merged.size() ; i++)
		assert(merged[i].get() == merged2[i].get());

	return 0;
}

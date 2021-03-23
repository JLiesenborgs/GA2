#include "population.h"
#include "mersennerandomnumbergenerator.h"
#include "vectorgenomefitness.h"
#include "basicnondominatedsetcreator.h"
#include <iostream>
#include <algorithm>

using namespace std;
using namespace errut;
using namespace mogal2;

int main(int argc, char const *argv[])
{
    random_device rd;
    MersenneRandomNumberGenerator rng(rd());

    Population pop;
    size_t N = 50;
    size_t dim = 2;
    
    for (size_t i = 0 ; i < N ; i++)
    {
        VectorGenome<float> dummy;
        VectorFitness<float> f { dim };
        f.setCalculated();

        for (size_t j = 0 ; j < dim ; j++)
            f.getValues()[j] = (int)(rng.getRandomDouble()*100);

        pop.append(make_shared<Individual>(dummy.createCopy(), f.createCopy()));
    }

    vector<shared_ptr<Individual>> ndSet;
    vector<shared_ptr<Individual>> remainder;

    BasicNonDominatedSetCreator bndsc(make_shared<VectorFitnessComparison<float>>(), dim);
    bool_t r;
    
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

    // if (!(r = bndsc.calculateNonDomitatedSet(pop.individuals(), ndSet, remainder)))
    //     cerr << "Error: " << r.getErrorString() << endl;

    printSet(pop.individuals());
    cout << endl << endl;
    // printSet(ndSet);
    // cout << endl << endl;
    // printSet(remainder);
    // cout << endl << endl;

    if (!(r = bndsc.calculateAllNDSets(pop.individuals())))
    {
        cerr << "Error: " << r.getErrorString() << endl;
        return -1;
    }

    cout << "# Found " << bndsc.getNumberOfSets() <<  " sets" << endl;
    for (size_t i = 0 ; i < bndsc.getNumberOfSets() ; i++)
    {
        auto s = bndsc.getSet(i);
        sort(s.begin(), s.end(), [](auto &e1, auto &e2) {
            const VectorFitness<float> &f1 = dynamic_cast<const VectorFitness<float> &>(e1->fitnessRef());
            const VectorFitness<float> &f2 = dynamic_cast<const VectorFitness<float> &>(e2->fitnessRef());
            return f1.getValues()[0] < f2.getValues()[0];
        });
        printSet(s);
        cout << endl << endl;
    }
    return 0;
}

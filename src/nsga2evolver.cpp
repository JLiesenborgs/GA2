#include "nsga2evolver.h"
#include "selection.h"
#include "tournamentparentselection.h"
#include "remainingtargetpopulationsizeiteration.h"
#include <numeric>
#include <algorithm>
#include <iostream>

#include "fasternondominatedsetcreator.h"

using namespace std;
using namespace errut;

namespace eatk
{

class DummySelPop : public SelectionPopulation
{
public:
	// Avoid bailing on any kind of check
	errut::bool_t check(const Population &population) { return true; }
	errut::bool_t processPopulation(const std::shared_ptr<Population> &population, size_t targetPopulationSize)
	{
		return true;
	}
};

NSGA2Evolver::NSGA2Evolver(
	const std::shared_ptr<RandomNumberGenerator> &rng,
	const std::shared_ptr<GenomeCrossover> &genomeCrossover,
	const std::shared_ptr<GenomeMutation> &genomeMutation,
	const std::shared_ptr<FitnessComparison> &fitComp, size_t numObjectives
	)
	: m_fitComp(fitComp), m_numObjectives(numObjectives)
{
	auto fitWrapComp = make_shared<FitWrapperComparison>();
	double cloneFrac = 0; // TODO: allow cloning?
	m_crossover = make_unique<SinglePopulationCrossover>(cloneFrac, true, 
		make_shared<DummySelPop>(),
		make_shared<TournamentParentSelection>(rng, 2, fitWrapComp),
		genomeCrossover, genomeMutation, nullptr,
		make_shared<RemainingTargetPopulationSizeIteration>(), rng);

	m_tmpPop = make_shared<Population>();
}

NSGA2Evolver::~NSGA2Evolver()
{

}

bool_t NSGA2Evolver::check(const std::shared_ptr<Population> &population)
{
	if (m_numObjectives < 2)
		return "Need at least two fitness objectives";

	if (population->size() < 2)
		return "Population size too small";

	for (auto &ind : population->individuals())
	{
		const Fitness &f = ind->fitnessRef();
		if (!f.hasRealValues())
			return "Fitness must consist of real values";
	}

	buildWrapperPopulation(*population);
	m_tmpPop->clear();
	for (auto &i : m_popWrapper)
		m_tmpPop->append(i);

	bool_t r;
	if (!(r = m_crossover->check(m_tmpPop)))
		return "Error checking final crossover: " + r.getErrorString();

	// TODO! Need to check anything else?

	return true;
}

void NSGA2Evolver::buildWrapperPopulation(const Population &population)
{
	m_popWrapper.clear();
	for (size_t i = 0 ; i < population.size() ; i++)
	{
		auto &ind = population.individual(i);
		m_popWrapper.push_back(make_shared<IndWrapper>(m_numObjectives, ind->genome(),
		                                               make_shared<FitWrapper>(ind->fitness()),
													   ind->getIntroducedInGeneration(), i));
	}
}

bool_t NSGA2Evolver::createNewPopulation(size_t generation, std::shared_ptr<Population> &population, size_t targetPopulationSize)
{
	// First generation (P_0):
	//   Need to extend population with Q_0, but can't do much more on this step
	//   since the fitness values need to be calculated
	// Next generation: we have R_0 = union of P_0 and Q_0, select P_1
	// (this already has fitness values) and generate Q_1 (doesn't have fitness
	// values yet)
	// Next generation, same as before

	// TODO: pass ndset creator as argument to constructor
	auto ndsCreator = make_shared<FasterNonDominatedSetCreator>(m_fitComp, m_numObjectives);
	bool_t r;

	// Create wrapper population, so that we can keep track of extra information
	buildWrapperPopulation(*population);

//	cout << "POPULATION[" << generation << "] = \n";
//	population->print();

	const Individual &refInd = *(population->individual(0));

	auto getIndividualFromWrapper = [&population,&refInd](const shared_ptr<Individual> &i1)
	{
		assert(dynamic_cast<const IndWrapper*>(i1.get()));
		const IndWrapper *pI1 = static_cast<const IndWrapper*>(i1.get());

		assert(pI1->m_originalPosition < population->size());
		return population->individual(pI1->m_originalPosition);
	};
	
	auto saveBestIndividuals = [this, &ndsCreator, &population, getIndividualFromWrapper](const vector<shared_ptr<Individual>> &ndSet)
	{
		const Individual &refInd = *(population->individual(0));
		m_best.clear();

		for (auto &ind : ndSet)
			m_best.push_back(getIndividualFromWrapper(ind)->createCopy());

		// TODO: remove duplicates necessary?
	};

	auto createNewPopulationFromNDRankAndCrowding = [this, targetPopulationSize, generation]() -> bool_t
	{
		assert(m_popWrapper.size() == targetPopulationSize);

		m_tmpPop->clear();
		for (auto &i : m_popWrapper)
			m_tmpPop->append(i);

		bool_t r;

		if (!(r = m_crossover->createNewPopulation(generation, m_tmpPop, targetPopulationSize)))
			return "Can't calculate genome crossover: " + r.getErrorString();

		m_popWrapper = m_tmpPop->individuals();

		return true;
	};

	auto unwrapPopulation = [this, targetPopulationSize, &population]()
	{
		assert(m_popWrapper.size() == targetPopulationSize*2);

		// Go from the wrappers back to the 'regular' population
		shared_ptr<Individual> refInd = population->individual(0); // Get one individual as reference
		population->clear();

		for (size_t i = 0 ; i < m_popWrapper.size() ; i++)
		{
			const shared_ptr<Individual> &wrapperInd = m_popWrapper[i];
			const shared_ptr<Fitness> &wrapperFit = wrapperInd->fitness();
			assert(dynamic_cast<FitWrapper*>(wrapperFit.get()));

			const FitWrapper *pWrapperFit = static_cast<const FitWrapper*>(wrapperFit.get());

			auto newInd = refInd->createNew(wrapperInd->genome(), pWrapperFit->m_origFitness, wrapperInd->getIntroducedInGeneration());
			newInd->setLastMutationGeneration(wrapperInd->getLastMutationGeneration());
#ifndef NDEBUG
			// First half should already have calculated fitnesses, last half shouldn't
			if (i < targetPopulationSize)
				assert(pWrapperFit->m_origFitness->isCalculated());
			else
				assert(!pWrapperFit->m_origFitness->isCalculated());
#endif
			population->append(newInd);
		}
	};

	// TODO: stop if targetPopulationSize has been reached/exceeded
	auto createNDSetsAndCalculateCrowdingDistances = [this,&ndsCreator]() -> bool_t
	{
		bool_t r;
		// Using wrapper to store original index positions
		if (!(r = ndsCreator->calculateAllNDSets(m_popWrapper)))
			return "Can't create non-dominated sets: " + r.getErrorString();
		
		// Store the ndset index for each individual
		size_t numSets = ndsCreator->getNumberOfSets();
		for (size_t s = 0 ; s < numSets ; s++)
		{
			for (auto &ind : ndsCreator->getSet(s))
			{
				assert(dynamic_cast<IndWrapper*>(ind.get()));
				IndWrapper *pIndWrapper = static_cast<IndWrapper*>(ind.get());

				assert(dynamic_cast<FitWrapper*>(pIndWrapper->fitnessPtr()));
				FitWrapper &fitness = static_cast<FitWrapper&>(pIndWrapper->fitnessRef());
				fitness.m_ndSetIndex = s;
			}

			// Also calculate crowding distances per set
			calculateCrowdingDistances(ndsCreator->getSet(s));
		}	
		return true;
	};

	auto getCrowdingValue = [](const shared_ptr<Individual> &i1) {
		assert(dynamic_cast<const IndWrapper*>(i1.get()));
		const IndWrapper *pI1 = static_cast<const IndWrapper*>(i1.get());
		const Fitness &f1 = pI1->fitnessRef();
		assert(dynamic_cast<const FitWrapper*>(&f1));
		const FitWrapper &fw = static_cast<const FitWrapper&>(f1);

		return fw.m_totalFitnessDistance;
	};

	if (population->size() == targetPopulationSize)
	{
		if (generation != 0)
			return "Expecting this population size only on generation 0";

		if (!(r = createNDSetsAndCalculateCrowdingDistances()))
			return r;

		saveBestIndividuals(ndsCreator->getSet(0));

		if (!(r = createNewPopulationFromNDRankAndCrowding()))
			return r;

		unwrapPopulation();
	}
	else if (population->size() == targetPopulationSize*2)
	{
		if (generation == 0)
			return "Unexpected double population size for generation 0";

		if (!(r = createNDSetsAndCalculateCrowdingDistances()))
			return r;

		//cerr << "Number of sets: " << ndsCreator->getNumberOfSets() << endl;

		vector<shared_ptr<Individual>> wrappersToKeep;
		// Check ND sets, copy until popsize filled, for last ND set
		// sort on crowing distance and keep fill the remaining pop size with
		// less crowded ones
		for (size_t setIdx = 0 ; setIdx < ndsCreator->getNumberOfSets() ; setIdx++)
		{
			auto &ndSetRef = ndsCreator->getSet(setIdx);
			if (wrappersToKeep.size() + ndSetRef.size() <= targetPopulationSize)
			{
				// Can keep entire set
				for (auto &i : ndSetRef)
					wrappersToKeep.push_back(i);

				if (setIdx == 0)
					saveBestIndividuals(ndSetRef);
			}
			else
			{
				vector<shared_ptr<Individual>> ndSetCopy = ndSetRef;

				// sort ndset on crowding distance, keep less crowded ones
				sort(ndSetCopy.begin(), ndSetCopy.end(), [&getCrowdingValue](const shared_ptr<Individual> &i1, const shared_ptr<Individual> &i2)
				{
					return getCrowdingValue(i1) > getCrowdingValue(i2);
				});

				ndSetCopy.resize(targetPopulationSize-wrappersToKeep.size());

				for (auto &i : ndSetCopy)
					wrappersToKeep.push_back(i);

				if (setIdx == 0)
					saveBestIndividuals(ndSetCopy);

				break; // Need to stop here
			}
		}

		swap(wrappersToKeep, m_popWrapper);

		if (!(r = createNewPopulationFromNDRankAndCrowding()))
			return r;

		unwrapPopulation();
	}
	else
		return "Unexpected population size " + to_string(population->size()) + " for target size " + to_string(targetPopulationSize);

	return true;
}

void NSGA2Evolver:: calculateCrowdingDistances(const std::vector<std::shared_ptr<Individual>> &ndset0) const
{
	auto ndset = ndset0;
	
	auto processObjective = [&ndset](size_t objectiveNumber)
	{
		auto getFitnessValue = [objectiveNumber](const shared_ptr<Individual> &i1) {
			assert(dynamic_cast<const IndWrapper*>(i1.get()));
			const IndWrapper *pI1 = static_cast<const IndWrapper*>(i1.get());
			const Fitness &f1 = pI1->fitnessRef();

			assert(f1.hasRealValues());
			return f1.getRealValue(objectiveNumber);
		};

		auto setDistanceValue = [objectiveNumber](shared_ptr<Individual> &i1, double value) {
			assert(dynamic_cast<IndWrapper*>(i1.get()));
			IndWrapper *pI1 = static_cast<IndWrapper*>(i1.get());

			assert(objectiveNumber < pI1->m_fitnessDistances.size());
			pI1->m_fitnessDistances[objectiveNumber] = value;
		};

		// Sort on fitness values for this objective
		sort(ndset.begin(), ndset.end(), [&getFitnessValue, objectiveNumber](const shared_ptr<Individual> &i1, const shared_ptr<Individual> &i2) {
			return getFitnessValue(i1) > getFitnessValue(i2);
		});

		assert(ndset.size() > 0);
		double maxFitness = getFitnessValue(ndset.front());
		double minFitness = getFitnessValue(ndset.back());
		double fDiff = maxFitness - minFitness;
		if (fDiff == 0)
			fDiff = 1.0; // avoid div by zero

		// Distance for this fitness for best and worst is infinity
		setDistanceValue(ndset.front(), numeric_limits<double>::infinity());
		setDistanceValue(ndset.back(), numeric_limits<double>::infinity());

		// Fill in the rest
		for (size_t i = 1 ; i < ndset.size()-1 ; i++)
		{
			double fNext = getFitnessValue(ndset[i+1]);
			double fPrev = getFitnessValue(ndset[i-1]);

			double dist = (fPrev - fNext)/fDiff;
			assert(dist >= 0);
			setDistanceValue(ndset[i], dist);
		}
	};

	for (size_t idx = 0 ; idx < m_numObjectives ; idx++)
		processObjective(idx);

	for (auto &ind : ndset)
	{
		assert(dynamic_cast<IndWrapper*>(ind.get()));
		IndWrapper *pIndWrapper = static_cast<IndWrapper*>(ind.get());

		assert(pIndWrapper->m_fitnessDistances.size() == m_numObjectives);
		
		assert(dynamic_cast<FitWrapper*>(pIndWrapper->fitnessPtr()));
		FitWrapper &fitness = static_cast<FitWrapper&>(pIndWrapper->fitnessRef());

		fitness.m_totalFitnessDistance = accumulate(pIndWrapper->m_fitnessDistances.begin(),
		                                            pIndWrapper->m_fitnessDistances.end(), 0);
	}
}

}
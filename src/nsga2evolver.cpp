#include "nsga2evolver.h"
#include "selection.h"
#include "tournamentparentselection.h"
#include "remainingtargetpopulationsizeiteration.h"
#include "fasternondominatedsetcreator.h"
#include "basicnondominatedsetcreator.h"
#include <numeric>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace errut;

namespace eatk
{

class DummySelPop : public SelectionPopulation
{
public:
	// Avoid bailing on any kind of check
	bool_t check(const Population &population) { return true; }
	bool_t processPopulation(const shared_ptr<Population> &population, size_t targetPopulationSize)
	{
		return true;
	}
};

NSGA2IndividualWrapper::NSGA2IndividualWrapper(size_t numObjectives,
		const shared_ptr<Genome> &genome, const shared_ptr<Fitness> &fitness,
		size_t introducedInGeneration, size_t originalPosition)
	: Individual(genome, fitness, introducedInGeneration),
		m_originalPosition(originalPosition)
{
	m_fitnessDistances.resize(numObjectives);
}

shared_ptr<Individual> NSGA2IndividualWrapper::createNew(const shared_ptr<Genome> &genome,
	const shared_ptr<Fitness> &fitness,
	size_t introducedInGeneration) const
{
	return make_shared<NSGA2IndividualWrapper>(m_fitnessDistances.size(), genome, fitness, introducedInGeneration, numeric_limits<size_t>::max());
}

string NSGA2IndividualWrapper::toString() const
{
	string s;

	s += "pos: " + to_string(m_originalPosition) + " |";
	for (auto d : m_fitnessDistances)
		s += " " + to_string(d);
	s += "| " + Individual::toString();
	return s;
}

NSGA2FitnessWrapper::NSGA2FitnessWrapper(const shared_ptr<Fitness> &origFitness)
	: m_origFitness(origFitness)
{
	m_totalFitnessDistance = numeric_limits<double>::quiet_NaN();
	m_ndSetIndex = numeric_limits<size_t>::max();
}

shared_ptr<Fitness> NSGA2FitnessWrapper::createCopy(bool copyContents) const
{
	shared_ptr<Fitness> realFitness = m_origFitness->createCopy(copyContents);
	shared_ptr<NSGA2FitnessWrapper> newFit = make_shared<NSGA2FitnessWrapper>(realFitness);
	if (copyContents)
	{
		newFit->m_totalFitnessDistance = m_totalFitnessDistance;
		newFit->m_ndSetIndex = m_ndSetIndex;
	}
	return newFit;
}

string NSGA2FitnessWrapper::toString() const
{
	return "{ ndset: " + to_string(m_ndSetIndex) + ", dist: " + to_string(m_totalFitnessDistance) + ", orig: " + m_origFitness->toString() + " }";
}

NSGA2FitnessWrapperOriginalComparison::NSGA2FitnessWrapperOriginalComparison(const shared_ptr<FitnessComparison> &origCmp)
	: m_origCmp(origCmp) 
{
}

bool_t NSGA2FitnessWrapperOriginalComparison::check(const Fitness &f) const
{
	if (!dynamic_cast<const NSGA2FitnessWrapper*>(&f))
		return "Expecting NSGA2FitnessWrapper object";
	const NSGA2FitnessWrapper &fw = static_cast<const NSGA2FitnessWrapper&>(f);
	return m_origCmp->check(*(fw.m_origFitness));
}

bool NSGA2FitnessWrapperOriginalComparison::isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const
{
	assert(dynamic_cast<const NSGA2FitnessWrapper*>(&first0) && dynamic_cast<const NSGA2FitnessWrapper*>(&second0));
	const NSGA2FitnessWrapper &first = static_cast<const NSGA2FitnessWrapper&>(first0);
	const NSGA2FitnessWrapper &second = static_cast<const NSGA2FitnessWrapper&>(second0);

	return m_origCmp->isFitterThan(*(first.m_origFitness), *(second.m_origFitness), objectiveNumber);
}

bool_t NSGA2FitWrapperNDSetCrowdingComparison::check(const Fitness &f) const
{
	if (!dynamic_cast<const NSGA2FitnessWrapper*>(&f))
		return "Expecting NSGA2FitnessWrapper object";
	return true;
}

bool NSGA2FitWrapperNDSetCrowdingComparison::isFitterThan(const Fitness &first0, const Fitness &second0, size_t objectiveNumber) const
{
	assert(dynamic_cast<const NSGA2FitnessWrapper*>(&first0) && dynamic_cast<const NSGA2FitnessWrapper*>(&second0));
	const NSGA2FitnessWrapper &first = static_cast<const NSGA2FitnessWrapper&>(first0);
	const NSGA2FitnessWrapper &second = static_cast<const NSGA2FitnessWrapper&>(second0);

	assert(first.m_ndSetIndex != numeric_limits<size_t>::max());
	assert(second.m_ndSetIndex != numeric_limits<size_t>::max());
	if (first.m_ndSetIndex < second.m_ndSetIndex)
		return true;
	if (second.m_ndSetIndex < first.m_ndSetIndex)
		return false;

	// Same ND set, check crowding distance
	assert(!isnan(first.m_totalFitnessDistance));
	assert(!isnan(second.m_totalFitnessDistance));

	// Use the less crowded one
	return first.m_totalFitnessDistance > second.m_totalFitnessDistance;
}

NSGA2Evolver::NSGA2Evolver(
	const std::shared_ptr<RandomNumberGenerator> &rng,
	const std::shared_ptr<GenomeCrossover> &genomeCrossover,
	const std::shared_ptr<GenomeMutation> &genomeMutation,
	const std::shared_ptr<FitnessComparison> &fitComp, size_t numObjectives
	)
	: m_numObjectives(numObjectives)
{
	m_fitOrigComp = make_shared<NSGA2FitnessWrapperOriginalComparison>(fitComp);
	auto fitWrapComp = make_shared<NSGA2FitWrapperNDSetCrowdingComparison>();
	double cloneFrac = 0; // TODO: allow cloning?
	m_crossover = make_unique<SinglePopulationCrossover>(cloneFrac, true, 
		make_shared<DummySelPop>(),
		make_shared<TournamentParentSelection>(rng, 2, fitWrapComp),
		genomeCrossover, genomeMutation, nullptr,
		make_shared<RemainingTargetPopulationSizeIteration>(), rng);

	m_tmpPop = make_shared<Population>();
	m_ndSetCreator = allocatedNDSetCreator(m_fitOrigComp, m_numObjectives);
}

NSGA2Evolver::~NSGA2Evolver()
{
}

shared_ptr<NonDominatedSetCreator> NSGA2Evolver::allocatedNDSetCreator(const std::shared_ptr<FitnessComparison> &fitCmp, size_t numObjectives)
{
	return make_shared<FasterNonDominatedSetCreator>(fitCmp, numObjectives);
	//return make_shared<BasicNonDominatedSetCreator>(fitCmp, numObjectives);
}

bool_t NSGA2Evolver::check(const std::shared_ptr<Population> &population)
{
	if (m_numObjectives < 2)
		return "Need at least two fitness objectives";

	if (population->size() < 2)
		return "Population size too small";

	if (!m_ndSetCreator)
		return "Need a non-dominated set creator!";

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

	return true;
}

void NSGA2Evolver::buildWrapperPopulation(const Population &population)
{
	m_popWrapper.clear();
	for (size_t i = 0 ; i < population.size() ; i++)
	{
		auto &ind = population.individual(i);
		m_popWrapper.push_back(make_shared<NSGA2IndividualWrapper>(m_numObjectives, ind->genome(),
		                                               make_shared<NSGA2FitnessWrapper>(ind->fitness()),
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

	assert(m_ndSetCreator.get());

	// Create wrapper population, so that we can keep track of extra information
	buildWrapperPopulation(*population);

//	cout << "POPULATION[" << generation << "] = \n";
//	population->print();

	const Individual &refInd = *(population->individual(0));

	auto getIndividualFromWrapper = [&population,&refInd](const shared_ptr<Individual> &i1)
	{
		assert(dynamic_cast<const NSGA2IndividualWrapper*>(i1.get()));
		const NSGA2IndividualWrapper *pI1 = static_cast<const NSGA2IndividualWrapper*>(i1.get());

		assert(pI1->m_originalPosition < population->size());
		return population->individual(pI1->m_originalPosition);
	};
	
	auto saveBestIndividuals = [this, &population, getIndividualFromWrapper](const vector<shared_ptr<Individual>> &ndSet)
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
			assert(dynamic_cast<NSGA2FitnessWrapper*>(wrapperFit.get()));

			const NSGA2FitnessWrapper *pWrapperFit = static_cast<const NSGA2FitnessWrapper*>(wrapperFit.get());

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

	auto createNDSetsAndCalculateCrowdingDistances = [this, targetPopulationSize]() -> bool_t
	{
		bool_t r;
		// Using wrapper to store original index positions
		if (!(r = m_ndSetCreator->calculateAllNDSets(m_popWrapper, targetPopulationSize)))
			return "Can't create non-dominated sets: " + r.getErrorString();
		
		// Store the ndset index for each individual
		size_t numSets = m_ndSetCreator->getNumberOfSets();
		for (size_t s = 0 ; s < numSets ; s++)
		{
			for (auto &ind : m_ndSetCreator->getSet(s))
			{
				assert(dynamic_cast<NSGA2IndividualWrapper*>(ind.get()));
				NSGA2IndividualWrapper *pIndWrapper = static_cast<NSGA2IndividualWrapper*>(ind.get());

				assert(dynamic_cast<NSGA2FitnessWrapper*>(pIndWrapper->fitnessPtr()));
				NSGA2FitnessWrapper &fitness = static_cast<NSGA2FitnessWrapper&>(pIndWrapper->fitnessRef());
				fitness.m_ndSetIndex = s;
			}

			// Also calculate crowding distances per set
			calculateCrowdingDistances(m_ndSetCreator->getSet(s));
		}	
		return true;
	};

	auto getCrowdingValue = [](const shared_ptr<Individual> &i1) {
		assert(dynamic_cast<const NSGA2IndividualWrapper*>(i1.get()));
		const NSGA2IndividualWrapper *pI1 = static_cast<const NSGA2IndividualWrapper*>(i1.get());
		const Fitness &f1 = pI1->fitnessRef();
		assert(dynamic_cast<const NSGA2FitnessWrapper*>(&f1));
		const NSGA2FitnessWrapper &fw = static_cast<const NSGA2FitnessWrapper&>(f1);

		return fw.m_totalFitnessDistance;
	};

	bool_t r;
	if (population->size() == targetPopulationSize)
	{
		if (generation != 0)
			return "Expecting this population size only on generation 0";

		if (!(r = createNDSetsAndCalculateCrowdingDistances()))
			return r;

		/*
		cerr << "Initial population" << endl;
		for (auto &i : m_popWrapper)
		{
			cerr << i->toString() << endl;
		}
		throw runtime_error("TODO");
		*/

		saveBestIndividuals(m_ndSetCreator->getSet(0));

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
		for (size_t setIdx = 0 ; setIdx < m_ndSetCreator->getNumberOfSets() ; setIdx++)
		{
			auto &ndSetRef = m_ndSetCreator->getSet(setIdx);
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
			assert(dynamic_cast<const NSGA2IndividualWrapper*>(i1.get()));
			const NSGA2IndividualWrapper *pI1 = static_cast<const NSGA2IndividualWrapper*>(i1.get());
			const Fitness &f1 = pI1->fitnessRef();

			assert(f1.hasRealValues());
			return f1.getRealValue(objectiveNumber);
		};

		auto setDistanceValue = [objectiveNumber](shared_ptr<Individual> &i1, double value) {
			assert(dynamic_cast<NSGA2IndividualWrapper*>(i1.get()));
			NSGA2IndividualWrapper *pI1 = static_cast<NSGA2IndividualWrapper*>(i1.get());

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
		assert(dynamic_cast<NSGA2IndividualWrapper*>(ind.get()));
		NSGA2IndividualWrapper *pIndWrapper = static_cast<NSGA2IndividualWrapper*>(ind.get());

		assert(pIndWrapper->m_fitnessDistances.size() == m_numObjectives);
		
		assert(dynamic_cast<NSGA2FitnessWrapper*>(pIndWrapper->fitnessPtr()));
		NSGA2FitnessWrapper &fitness = static_cast<NSGA2FitnessWrapper&>(pIndWrapper->fitnessRef());

		// Note: this has to be 0.0 and not just 0, otherwise the doubles will be cast to
		//       int values
		fitness.m_totalFitnessDistance = accumulate(pIndWrapper->m_fitnessDistances.begin(),
		                                            pIndWrapper->m_fitnessDistances.end(), 0.0);
	}
}

}
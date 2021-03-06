//#include "krktestinstance.h"
//#include "dataset.h"
#include "population.h"
#include <limits>

using namespace std;

// GLOBAL VARIABLES
// . . . .
//template <class T>
Population::Population(int n, int nprocs, DataSet<KRKTestInstance> ts, int popsize)
// n is the rank within the cluster of this island,
// and nprocs in the total number of compute-nodes in this run
// ts is the DataSet to use for evaluating this population's fitness
// popsize is how large the population will be
{
	myDataSet = ts;
	POP_SIZE = popsize;

	//printf("Population: starting constructor\n");
	//set the default values
	generation = 0;
	maxFitness = 0.0;

	b1 = new Individual2[POP_SIZE];  // need to randomize this one
	b2 = new Individual2[POP_SIZE];
	for (int i=0; i < POP_SIZE ; i++) b1[i].setRandomRule();
	mypop = b1;
	newpop = b2;
	newpop_count = 0;

	// set myrank and that of the neighbors for this island
	myrank = n;
	numprocs = nprocs;
	//printf("Population: finished constructor\n");
}

//template <class T>
Individual2* Population::getIndividual(int index) {
	return &mypop[index];
}

// computes the hamming distance between the input individual and each individual in this population
//template <class T>
vector<int> Population::calculateHammingForAll(Individual2* input){
	unsigned int * inputRule = (*input).getIntRule();
	unsigned int * compareRule;
	vector<int> hamming (POP_SIZE);
	for (int i = 0; i < POP_SIZE; ++i) {
		compareRule = mypop[i].getIntRule();
		hamming[i] = calculateHammingPair(inputRule, compareRule);
	}
	return hamming;
}

//  computes the hamming diversity within this population
//template <class T>
vector<float> Population::getInternalHammingDiversity(){
	unsigned int * iRule;
	unsigned int * jRule;
	int best = 0;
	int worst = 1000;
	int total = 0;
	vector< vector<int> > diversetable (POP_SIZE, vector<int>(POP_SIZE));
	// Get the best, worst, and total hamming distances. Also get the denominator for calculating the average.
	for (int i = 0; i < POP_SIZE; ++i) {
	  iRule = mypop[i].getIntRule();
	  for (int j = i+1; j < POP_SIZE; ++j) {
	    jRule = mypop[j].getIntRule();
	    int temp = calculateHammingPair(iRule, jRule);
	    diversetable[i][j] = temp;
	    diversetable[j][i] = temp;
	    best = max(best, temp);
	    worst = min(worst, temp); // In the case of internal diversity, we allow for worst to be 0 so that we can see uniformity among individuals.
	  }
	}

	//caclulate the total diversities for each individual
	vector<int> individs (POP_SIZE);
	int counter = 0;
	int ptotal = 0;
	for (int i=0; i < POP_SIZE; i++){
	  for (int j=0; j < POP_SIZE; j++){
	    individs[i] += diversetable[i][j];
	    ptotal += diversetable[i][j];
	  }
	  individs[i] /= POP_SIZE - 1; // note: int division
	}

	//calculate the variance
	float variance = 0.0;
	float average = (float) ptotal / (float) (POP_SIZE * (POP_SIZE - 1));
	for (int i = 0; i < POP_SIZE ; ++i) {
		variance += powf(((float) individs[i]/POP_SIZE - average), 2.0);
	}
	vector<float> output (4);
	output[0] = (float) best;
	output[1] = (float) worst;
	output[2] = average;
	output[3] = variance / (float) POP_SIZE;
	return output;
}

// Calculates the hamming distance for a pair of strings
//template <class T>
int Population::calculateHammingPair(unsigned int * a, unsigned int * b) {
	int HamDiff = 0;
	for (int i = 0; i < RULE_LEN; ++i) {
		if (a[i] != b[i]) ++HamDiff;
	}
	return HamDiff;
}

//template <class T>
void Population::addPopulationBitTotal(float * totals) {
	unsigned int * tempRule;
	for (int i = 0; i < POP_SIZE; ++i) {
		tempRule = mypop[i].getIntRule();
		for (int j = 0; j < RULE_LEN; ++j) {
			totals[j] += tempRule[j];
		}
	}
}

//template <class T>
void Population::updatePopulationIntRules() {
	for (int i = 0; i < POP_SIZE; ++i) {
		mypop[i].resetIntRule();
	}
}

// update fitness of all individuals

// determine fitness with respect to the training data for breeding/survival/migration purposes
//template <class T>
void Population::updatePopulationFitness(char fitnessClassifier, char whichSet) {
	//cout << "update pop Fittenss:" << fitnessClassifier << endl;
	if (fitnessClassifier == ' ' || whichSet == ' ') exit(0);
	maxFitness = 0.0;
	totalFitness = 0.0;
	totalInverseFitness = 0.0;
	stdev = 0.0;
	bestIndiv = mypop[0];
	for (int i = 0; i < POP_SIZE; i++)
	{
		// printf("updating individual %d\n", i);
		updateFitness(&myDataSet, &mypop[i], fitnessClassifier, whichSet);
		//printf("updated individual %d\n", i);
		totalFitness += mypop[i].getFitness();
		totalInverseFitness += 1 - mypop[i].getFitness();
		if(maxFitness < mypop[i].getFitness())
		{
			maxFitness = mypop[i].getFitness();
			bestIndiv = mypop[i];
		}
	}

	avgFitness = totalFitness/POP_SIZE;
	for (int j = 0; j < POP_SIZE; j++)
	{
		float diff = mypop[j].getFitness() - avgFitness;
		stdev += diff * diff/(POP_SIZE - 1);
	}
	stdev = sqrt(stdev);
}

//template <class T>
void Population::updateFitness(DataSet<KRKTestInstance>* ts, Individual2* individual, char fitnessClassifier, char whichSet)
  /** Update the fitness of this individual over the set of krktestinstances using the appropriate fitness measure.
   */
{
  //printf("Entering updatFitness\n");
  float sum = 0.0;
  float testsofthistype;
  individual->setFitness(0.0);
  individual->resetConfMat();
  int setSize = ts->trainSetSize();
  if (whichSet == 'f') setSize = ts->testSetSize();
  KRKTestInstance krk;
  for(int i=0; i < setSize; i++) {
	  if (whichSet == 'f') krk = (KRKTestInstance)*ts->getTestI(i); // for full test set
	  else if (whichSet == 't') krk = (KRKTestInstance)*ts->getTrainI(i); // for training set
    switch (fitnessClassifier)
      {
      case 'h': {
    	  sum += krk.fitnessHiFi(individual); //cout << sum << " exiting" << endl; exit(0);
    	  break; }		// HiFi
      case 'l': {
    	  if(krk.classify(individual)) sum++;
    	  break; }// LoFi
      }
    //printf("On test case number %d\n", i);
    //printf("TC %d: inst: %s, classified as: %d, tc-binrep: %s\n", i, testSet[i].humanReadable().c_str(), classify(testSet[i]), testSet[i].getStringRep().c_str());
  }
  // alternate fitness computation based on the confusion matrix
  for(int i=0; i < RULE_CASES; i++)
    {
      testsofthistype = 0;
      for (int j=0; j < RULE_CASES; j++)
	testsofthistype += individual->confMat[i][j];
      individual->fitness += ((testsofthistype == 0) ? 0 : individual->confMat[i][i]/testsofthistype);
    }
  individual->fitness = individual->fitness/RULE_CASES;

  //fitness = sum / ts.NUM_TEST_CASES_TO_USE;
  //printf("Leaving updateFitness\n");
}

//template <class T>
void Population::updatePopulationRelevance(vector<float>* relevance) {
	for (int i = 0; i < POP_SIZE; ++i) {
		mypop[i].updateDiversityRelevance((*relevance)[i]);
	}
}

// select individuals for . . . .
//template <class T>
void Population::selectToSurvive(int n)
// Select individuals to survive to next population
{
	static Individual2 tmpI;
	int selectedIndividual;
	int availablepop = POP_SIZE;
	for (int i=0; i < n; i++)
	{
		//printf("Node %d: selecting individual to survive %d\n", myrank, i);
		selectedIndividual = selectIndividual(availablepop);
		//printf("Node %d: selected individual %d\n", myrank, selectedIndividual);
		if( !mypop[selectedIndividual].isSelected() )
		{
			newpop[newpop_count++] = mypop[selectedIndividual];
			mypop[selectedIndividual].select();
			tmpI = mypop[availablepop-1];
			mypop[--availablepop] = mypop[selectedIndividual];
			mypop[selectedIndividual] = tmpI;
		} else {
			printf("Node %d: FAILED TO SELECT AN UNSELECTED INDIVIDUAL! (try %d of %d)\n", myrank, i, n); exit(3);
		}
	}
	unselectAll();
}

//template <class T>
void Population::selectRandToMigrate(Individual2 * migrants, int num_migrants)
// stuffs randomly selected individuals (without replacement) into the migrants array
// ??? and sets their selected flag in the base population
{
	//cout << "\nRandom migration" << endl;
	int i;
	for(int n = 0; n < num_migrants; n++)
	{
		i = rand()%POP_SIZE;
		while(mypop[i].isSelected())  i=(i+1)%POP_SIZE;
		migrants[n] = mypop[i];
		//mypop[i].select();
	}
}

//template <class T>
void Population::selectStrongToMigrate(Individual2 * migrants, int num_migrants)
// stuffs individuals (bias towards strong) into the migrants array
// ??? and sets their selected flag in the base population
{
	int n =0;
	while(n < num_migrants)
	{
		int selectedIndividual = selectIndividual(POP_SIZE-n);
		if( !mypop[selectedIndividual].isSelected() )
		{
			//mypop[selectedIndividual].select();
			migrants[n] = mypop[selectedIndividual];
			n++;
		} else {
			printf("selectStrongToMigrate: selectIndividual returned an already selected individual %d\n", selectedIndividual);
		}
	}
}

//template <class T>
void Population::selectWeakToMigrate(Individual2 * migrants, int num_migrants)
// stuffs individuals (bias towards weak) into the migrants array and sets
// their selected flag in the base population
{
	int n =0;

	while(n < num_migrants)
	{
		int selectedIndividual = selectWeakIndividual();
		if( !mypop[selectedIndividual].isSelected() )
		{
			mypop[selectedIndividual].select();
			migrants[n] = mypop[selectedIndividual];
			n++;
		}
	}
}

// add immigrating individuals
//template <class T>
void Population::processImmigrants(Individual2 * v, int n)
{
	for (int i=0; i < n; i++)
		processOneImmigrant(v[i]);
}

//template <class T>
void Population::processOneImmigrant(Individual2 i)
{
	newpop[newpop_count++] = i;
}

// generate offspring
//template <class T>
void Population::generateOffspring(int n)
//breed remaining individuals two at a time
{
	//printf("Entering generateOffspring\n");
	//Individual2 * kids;
	Individual2 kids[2];
	for(int i=0; i < n; i+=2)
	{
		int parent1 = selectIndividual(POP_SIZE);
		int parent2 = selectIndividual(POP_SIZE);
		while( parent1 == parent2 ) {parent2 = selectIndividual(POP_SIZE); }
		//printf("Node %d(generateOffspring for %d of %d): selected two parents -- ready to breed(%d+%d)\n", myrank, i, n,  parent1, parent2);
		mypop[parent1].breedNCross(kids, mypop[parent2]);
		//printf("Node %d(generateOffspring): finished breeding\n", myrank);
		newpop[newpop_count++] = kids[0];
		//printf("added first kid\n");
		if(newpop_count < POP_SIZE) newpop[newpop_count++] = kids[1];
		//printf("added second kid, about do delete\n");
	}
	//printf("Leaving generateOffspring\n");
}

// switch to next generation
//template <class T>
void Population::nextGeneration(int PROB_MUTATE)
{
	//printf("Entering nextGeneration\n");
	Individual2 * tptr;
	if (newpop_count != POP_SIZE) {printf("DID NOT FILL NEWPOPULATION\n"); exit(-2);}

	// mutate maybe
	for(int i=0; i < POP_SIZE; i++)
		if (rand() < PROB_MUTATE)
			newpop[i].mutate();

	// move newpop to mypop
	// delete [] mypop;
	tptr = mypop;
	mypop = newpop;
	//newpop = new Individual2[POP_SIZE]; // for next generation
	newpop = tptr;
	newpop_count = 0;

	generation++;

	// unselect all individuals
	//unselectAll(); ** this is done after selectToSurvive and no other selection setting taking place
	//printf("Leaving nextGeneration\n");
	return;
}

// private member functions
//template <class T>
void Population::unselectAll()
{
	Individual2 * tpop_ptr = mypop;
	for(int j = 0; j < POP_SIZE; j++)
	{
		tpop_ptr->unselect();
		tpop_ptr++;
	}
}

//template <class T>
int Population::selectIndividual(int availablepop)
{
	switch (WHICH_SELECT) {
	case 0: return rand() % availablepop; break;// completely random selection with no bias either toward fitness or diversity
	case 1: return tournamentSelect(availablepop); break;
	case 2: return altSelectIndividual(availablepop); break;
	case 3: return relevanceTournamentSelect(availablepop); break;
	case 4: return relevanceTournamentSelect(availablepop); break;
	}
}

//template <class T>
int Population::relevanceTournamentSelect(int availablepop) {
	int bestIndex, candidate;
	bestIndex = rand() % availablepop;
	float bestFit = mypop[bestIndex].getDiversityRelevance();
	//cout << " BI: " << bestIndex << " V: " << bestFit << endl;
	for (int i=1; i < TOURNAMENT_SIZE; i++) {
		candidate = rand() % availablepop;
		if (mypop[candidate].getDiversityRelevance() > bestFit) {
			bestIndex = candidate;
			bestFit = mypop[candidate].getDiversityRelevance();
			//cout << " BI: " << bestIndex << " V: " << bestFit << endl;
		}
	}
	return bestIndex;
}

//template <class T>
int Population::tournamentSelect(int availablepop)
// tournament selection with replacement for tournament participants -- best of tournament candidates selected
{
	float bestFit;
	int bestIndex, candidate;
	bestIndex = rand() % availablepop;
	bestFit = mypop[bestIndex].getFitness();
	for (int i=1; i < TOURNAMENT_SIZE; i++) {
		candidate = rand() % availablepop;
		if (mypop[candidate].getFitness() > bestFit) {
			bestIndex = candidate;
			bestFit = mypop[candidate].getFitness();
		}
	}
	return bestIndex;
}

//template <class T>
int Population::altSelectIndividual(int availablepop)
{
	float epsilon = std::numeric_limits<float>::epsilon();
	static float MY_RAND_MAX = (float)RAND_MAX + 1.0;
	float availableFitness = 0.0;
	//printf("Node %d(altSelectIndividual): About to accumulate available fitness\n", myrank);
	for(int i=0; i < availablepop; i++)
	  availableFitness += mypop[i].getFitness();
	//printf("Node %d(altSelectIndividual): Finished accumulating available fitness\n", myrank);
	//float arandnum = ((float)RAND_MAX/MY_RAND_MAX) * (totalFitness - usedFit); // for testing floating point precision problem
	float arandnum = (((float)rand())/MY_RAND_MAX) * availableFitness;
	//arandnum -= (arandnum * epsilon); I don't think this is necessary
	int selectedIndex = 0;
	while ( (selectedIndex < availablepop-1) && (arandnum > mypop[selectedIndex].getFitness() ))
	{
	  arandnum -= mypop[selectedIndex].getFitness();
	  selectedIndex++;
	}
	//if ( selectedIndex >= availablepop ) {
	//   printf("altSelectIndividual: COULD NOT SELECT INDIVIDUAL (MY_RAND_MAX %f; arandnum remaining %f; selectedIndex %i;\n", MY_RAND_MAX, arandnum, selectedIndex); exit(-3);}
	//cout << selectedIndex << endl;
	return selectedIndex;
}

//template <class T>
int Population::selectWeakIndividual()
{
	float arandnum = (((float)rand())/RAND_MAX) * totalInverseFitness;
	int selectedIndex = 0;
	Individual2 * iptr = mypop;
	while ( (arandnum > (1 - iptr->getFitness())) && (selectedIndex < POP_SIZE) )
	{
		arandnum -= (1 - iptr->getFitness());
		iptr++;
		selectedIndex++;
	}
	return selectedIndex;
}

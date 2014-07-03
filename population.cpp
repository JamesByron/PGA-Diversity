#include "testinstance2.h"
#include "testset.h"
#include "individual2.h"
#include "population.h"


using namespace std;

// GLOBAL VARIABLES
// . . . .


Population::Population(int n, int nprocs, TestSet ts, int popsize)
  // n is the rank within the cluster of this island,
  // and nprocs in the total number of compute-nodes in this run
  // ts is the TestSet to use for evaluating this population's fitness
  // popsize is how large the population will be
{
  myTestSet = ts;
  POP_SIZE = popsize;
  HEMMING_DIST = new int[POP_SIZE][POP_SIZE];

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

int calculateHemmingPair(string a, string b) {
	int HemDiff = 0;
	if (a.length() != b.length) return -1;	// complain
	for (int i = 0; i < a.length(); ++i) {
		if (a[i] != b[i]) ++HemDiff;
	}
	return HemDiff;
}

void Population::updateHemmingDist() {
	for (int i = 0; i < POP_SIZE; ++i) {
		string iString = mypop[i].getStringRule();
		for (int j = i+1; j < POP_SIZE; ++j) {
			string jString = mypop[j].getStringRule();
			HEMMING_DIST[i][j] = calculateHemmingPair(iString, jString);
		}
	}
}


// update fitness of all individuals

// determine fitness with respect to the training data for breeding/survival/migration purposes
void Population::updatePopulationFitness(char WHICH_FITNESS)
{
  maxFitness = 0.0;
  totalFitness = 0.0;
  totalInverseFitness = 0.0;
  stdev = 0.0;
  bestIndiv = mypop[0];
  for (int i = 0; i < POP_SIZE; i++)
  {
    // printf("updating individual %d\n", i);
    mypop[i].updateFitness(myTestSet, WHICH_FITNESS);
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

void Population::updatePopulationFitness(vector<TestInstance2> allti, char WHICH_CLASSIFY)
{
  maxFitness = 0.0;
  totalFitness = 0.0;
  totalInverseFitness = 0.0;
  stdev = 0.0;
  bestIndiv = mypop[0];
  for (int i = 0; i < POP_SIZE; i++)
  {
    mypop[i].updateFitness(allti, WHICH_CLASSIFY);
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

// determine fitness as measure of accuracy over the actual holdout test set
void Population::populationAccuracy(char WHICH_CLASSIFY)
{
  maxFitness = 0.0;
  totalFitness = 0.0;
  stdev = 0.0;
  bestIndiv = mypop[0];
  for (int i = 0; i < POP_SIZE; i++)
  {
    mypop[i].findAccuracy(myTestSet, WHICH_CLASSIFY);
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

// select individuals for . . . .

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

void Population::selectRandToMigrate(Individual2 * migrants, int num_migrants)
  // stuffs randomly selected individuals (without replacement) into the migrants array 
  // ??? and sets their selected flag in the base population
{
  int i;
  for(int n = 0; n < num_migrants; n++)
  {
    i = rand()%POP_SIZE;
    while(mypop[i].isSelected())  i=(i+1)%POP_SIZE;
    migrants[n] = mypop[i];
    //mypop[i].select();
  }
}

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

void Population::processImmigrants(Individual2 * v, int n)
{
  for (int i=0; i < n; i++)
    processOneImmigrant(v[i]);
}

void Population::processOneImmigrant(Individual2 i)
{
  newpop[newpop_count++] = i;
}

// generate offspring

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
      while( parent1 == parent2 ) { parent2 = selectIndividual(POP_SIZE); }
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


void Population::unselectAll()
{
  Individual2 * tpop_ptr = mypop;
  for(int j = 0; j < POP_SIZE; j++)
  {
    tpop_ptr->unselect();
    tpop_ptr++;
  }
}


int Population::selectIndividual(int availablepop)
{
  switch (WHICH_SELECT) {
  case 1: return tournamentSelect(availablepop); break;
  case 2: return altSelectIndividual(); break;
  }
}


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

int Population::altSelectIndividual()
{
  static float MY_RAND_MAX = (float)RAND_MAX + 1.0;
  float usedFit = 0.0;
  int usedCount = 0;
  //printf("Node %d(altSelectIndividual): About to accumulate used fitness\n", myrank);
  for(int i=0; i < POP_SIZE; i++)
    if ( mypop[i].isSelected() ) {
      usedFit += mypop[i].getFitness();
      usedCount++;
    }
  //printf("Node %d(altSelectIndividual): Finished accumulating used fitness\n", myrank);
  if ( usedCount == POP_SIZE ) { printf("Node %d(altSelectIndividual): ENTIRE POPULATION IS SELECTED!\n", myrank); exit(-3); }
  float arandnum = (((float)rand())/MY_RAND_MAX) * (totalFitness - usedFit);
  int selectedIndex = 0;
  while ( (selectedIndex < POP_SIZE) && (arandnum > mypop[selectedIndex].getFitness() || mypop[selectedIndex].isSelected() ) )
  {
    if ( !mypop[selectedIndex].isSelected() )
      arandnum -= mypop[selectedIndex].getFitness();
    selectedIndex++;
  }
  if ( selectedIndex >= POP_SIZE ) {printf("COULD NOT SELECT INDIVIDUAL\n"); exit(-3);}
  return selectedIndex;
}

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

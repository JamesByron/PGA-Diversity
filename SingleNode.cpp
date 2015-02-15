#include "testinstance2.h"
#include "testset.h"
#include "individual2.h"
#include "population.h"
#include "SingleNode.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

// GLOBAL VARIABLES
FILE *logFile1; // for Generations data
FILE *logFile2; // for Islands data
// migration
char WHICH_MIGRATION; // r: Random, s: Strong, w: Weak, n: None, f: Fake random migrant, . . .
int WHEN_MIGRATE;
int NUM_MIGRANTS_PER_ISLAND;
int NUM_IMMIGRANTS;
// topography
int NUM_ISLANDS;
int NUM_NEIGHBORS;
// population
int POP_SIZE;
int NUM_TEST_CASES_TO_USE;
int DEPTH_OF_TEST_INSTANCES;
// other
int MAX_GENERATIONS;
char WHICH_FITNESS;
char WHICH_CLASSIFY;
int PROB_MUTATE; // an int < RAND_MAX representing probability a single individual will get "hit"
int WHEN_FULL_TEST = 100;
// use RULE_LEN instead: int genomeLength = NUM_FEATURES*RULE_CASES*8;

SingleNode * islands;

SingleNode::SingleNode(int r, TestSet ts)
{
	myrank = r;
	//  printf("SingleNode: construct an instance of a single-node\n");
	//set the default values

	//  printf("SingleNode: about to create a population\n");
	a_pop = new Population(r, NUM_ISLANDS, ts, POP_SIZE);
	a_pop->updatePopulationIntRules();
	//  printf("SingleNode: have population and now update that pop's fitness\n");
	a_pop->updatePopulationFitness(WHICH_FITNESS);

	customs = new Individual2[NUM_IMMIGRANTS];
}

SingleNode::SingleNode()
// Default constructor with empty test-set
{

}

void SingleNode::doOneGeneration(int thisgen)
// do the stuff for one generation
{
	//printf("Node %d: sendMigrantStrings()\n", myrank);
	if (WHICH_MIGRATION != 'n' &&  !(thisgen % WHEN_MIGRATE) )
		sendMigrants();

	// survive
	//printf("Node %d: a_pop->selectToSurvive(%d)\n", myrank, ((int)(POP_SIZE*PROB_REMAIN)));
	a_pop->selectToSurvive((int)(POP_SIZE*PROB_REMAIN));

	// receive
	//printf("Node %d: receiveMigrantStrings()\n", myrank);
	// ignore when only one island, but need to receive migrants (unless sendMigrants takes care of it)
	if (WHICH_MIGRATION != 'n' && !(thisgen % WHEN_MIGRATE) )
		a_pop->processImmigrants(islands[myrank].customs, NUM_IMMIGRANTS);

	// breed
	//printf("Node %d: a_pop->generateOffspring(%d)\n", myrank,(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + NUM_IMMIGRANTS)));
	a_pop->generateOffspring(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + ((WHICH_MIGRATION != 'n' && !(thisgen % WHEN_MIGRATE) ) ? NUM_IMMIGRANTS : 0)));
	//a_pop->generateOffspring(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + NUM_IMMIGRANTS));

	// switch to next gen
	//printf("Node %d: a_pop->nextGeneration()\n", myrank);
	a_pop->nextGeneration(PROB_MUTATE);
	a_pop->updatePopulationIntRules();
	a_pop->updatePopulationFitness(WHICH_FITNESS);
}

void SingleNode::updateNodeIntRules() {
	a_pop->updatePopulationIntRules();
}

void SingleNode::updateNodeRelevance(vector<float>* islandRelevance) {
	a_pop->updatePopulationRelevance(islandRelevance);
}

Individual2* SingleNode::getIndividual(int index) {
	return a_pop->getIndividual(index);
}

void SingleNode::addIslandBitTotal(float * totals) {
	a_pop->addPopulationBitTotal(totals);
}

int SingleNode::sendMigrants()
{
	Individual2 * migrants = new Individual2[NUM_IMMIGRANTS];
	switch (WHICH_MIGRATION)
	{
	case 'r': a_pop->selectRandToMigrate(migrants, NUM_IMMIGRANTS); break;
	case 's': a_pop->selectStrongToMigrate(migrants, NUM_IMMIGRANTS); break;
	case 'w': a_pop->selectWeakToMigrate(migrants, NUM_IMMIGRANTS); break;
	case 'f': // fake migrants -- generate random individuals on receiving end
	case 'n': break;
	}

	for(int i=0; i < NUM_IMMIGRANTS; i++)
	{
		islands[(myrank+(1+(i%NUM_NEIGHBORS)))%NUM_ISLANDS].customs[i] = migrants[i];
	}

	delete [] migrants;
	return NUM_IMMIGRANTS;
}

// compares every individual to every other individual in all the nodes.
// This operation is very expensive.
vector<float> getPairwiseHammingDiversityForAllNodes() {
	int TOTAL_POP = POP_SIZE * NUM_ISLANDS;
	int bestDiversity = 0;
	int worstDiversity = 1000; // Start with a high number.
	vector<int> tempV;
	Individual2* currentIndividual;
	SingleNode currentNode;
	vector< vector<int> > diversetable (TOTAL_POP, vector<int>(TOTAL_POP));
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		currentNode = islands[i];
		for (int j = 0; j < POP_SIZE; ++j) {
			currentIndividual = currentNode.getIndividual(j);
			for (int k = i; k < NUM_ISLANDS; ++k) {
				tempV = islands[k].a_pop->calculateHammingForAll(currentIndividual);
				// now add these to the table
				for (int x=0; x < POP_SIZE; x++){
					diversetable[i*POP_SIZE+j][k*POP_SIZE+x] = tempV[x];
					diversetable[k*POP_SIZE+x][i*POP_SIZE+j] = tempV[x];
				}
			}
		}
	}

	// go back and compute each individual's average hamming distance distance
	vector<int> individs (POP_SIZE * NUM_ISLANDS);
	int ptotal = 0;
	for (int i=0; i < POP_SIZE*NUM_ISLANDS; i++){
		for (int j=0; j < POP_SIZE*NUM_ISLANDS; j++){
			individs[i] += diversetable[i][j];
			ptotal += diversetable[i][j];
			//if ((0 == diversetable[i][j]) && (i != j)) cout << "i:" << i << " j:" << j << endl; // Which of the individuals are identical
			if (bestDiversity < diversetable[i][j]) bestDiversity  = diversetable[i][j];
			else if ((worstDiversity > diversetable[i][j]) && (i != j)) worstDiversity = diversetable[i][j];
		}
		individs[i] /= POP_SIZE*NUM_ISLANDS - 1;
	}

	//calculate the variance
	float varianceTotal = 0.0;
	float average = (float) ptotal / (float) (POP_SIZE*NUM_ISLANDS * (POP_SIZE*NUM_ISLANDS - 1));
	for (int i = 0; i < POP_SIZE*NUM_ISLANDS; ++i) {
		varianceTotal += powf(((float) individs[i]/(POP_SIZE*NUM_ISLANDS) - average), 2.0);
	}
	vector<float> output (4);
	output[0] = (float) bestDiversity;
	output[1] = (float) worstDiversity;
	output[2] = average;
	output[3] = varianceTotal / ((float) POP_SIZE*NUM_ISLANDS);
	return output;
}

// for debugging purposes
void seeInsideVector(vector<float> input, string leadStr) {
	bool bad = false;
	cout << endl << leadStr;
	for (unsigned int i = 0; i < input.size(); ++i) {
		if ((i % 10) == 0) cout << endl;
		if (isnanf(input[i]) != 0) bad = true;
		cout << input[i] << ", ";
	}
	cout << endl;
	cout.flush();
	if (bad) {cout << endl << endl << "Found nan!" << endl; cout.flush(); exit(-1);}
}

void computeHammingBitTotals(float * bitTotals) {
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		islands[i].addIslandBitTotal(bitTotals);
	}
	int totalIndividuals = NUM_ISLANDS*POP_SIZE;
	for (int i = 0; i < RULE_LEN; ++i) {
		bitTotals[i] /= (float)totalIndividuals;
		//if (bitTotals[i] > 1.0) exit(-1);
	}
}

vector<float> getHammingRelevanceDetail() {
	float bitTotals[NUM_FEATURES*RULE_CASES*8] = {0.0};
	computeHammingBitTotals(bitTotals);
	int currentIndividual;
	vector<float> diversityStore (NUM_ISLANDS*POP_SIZE);
	unsigned int * genome;
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		for (int j = 0; j < POP_SIZE; ++j) {
			currentIndividual = i*POP_SIZE + j;
			diversityStore[currentIndividual] = 0;
			genome = (*islands[i].getIndividual(j)).getIntRule();
			for (int k = 0; k < RULE_LEN; ++k) {
				diversityStore[currentIndividual] += abs((float)genome[k] - bitTotals[k]);
			}
		}
	}
	return diversityStore;
}

void combineRelevanceWithFitness(float weight, vector<float>* relevance) {
	float inverseW = 1-weight;
	if ((weight > 1.0) || (weight < 0.0)) {cout << "combineRelevanceWithFitness: Weight must be between 0.0 and 1.0" << endl; exit(-1);}
	//if ((*relevance).size() != (NUM_ISLANDS*POP_SIZE)) {cout << "combineRelevanceWithFitness: relevance vector is of an unexpected size." << endl; exit(-1);}
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		for (int j = 0; j < POP_SIZE; ++j) {
			(*relevance)[(i*POP_SIZE)+j] = ((*relevance)[(i*POP_SIZE)+j]*weight) + ((*islands[i].getIndividual(j)).getFitness()*inverseW);
		}
	}
}

/*
 * Computes the variance of the fitness of the individuals in the entire population with the given test set.
 * Returned vector is <best, worst, average, variance>
 */
vector<float> getFitnessDiversity(int startIsland, int endIsland) {
	if ((startIsland < 0) || (endIsland > NUM_ISLANDS)) { cout << "Error in getFitnessDiversity; island index out of range " << startIsland << ", " << endIsland << endl; exit(-1); }
	int counter = 0;
	float worst = 1000.0;
	float mean, x, delta, var, best = 0.0;
	for (int i = startIsland; i < endIsland; ++i) {
		//islands[i].a_pop->updatePopulationFitness(testInstances, classification);
		for (int j = 0; j < POP_SIZE; ++j) {
			++counter;
			x = (*islands[i].getIndividual(j)).getFitness();
			if (x < worst) worst = x;
			else if (x > best) best = x;
			delta = x - mean;
			mean = mean + (delta / (float)counter);
			var = var + (delta * (x - mean));
		}
	}
	vector<float> output (4);
	output[0] = best;
	output[1] = worst;
	output[2] = mean; //the accuracy of mean is probably sufficient
	output[3] = var / (float)counter; // counter = POP_SIZE * NUM_ISLANDS
	return output;
}

/*
 * Computes the max, min, average, and variance of the numbers in the given vector.
 * Returned vector is <best, worst, average, variance, number of test instances that have 0 weight>
 */
vector<float> getDiversityValues(vector<float>* values, int numZeros) {
	int counter = 0;
	float worst = 99999999;
	float mean, x, delta, var, best = 0.0;
	for (int i = 0; i < (*values).size(); ++i) {
		++counter;
		x = (*values)[i];
		if (x < worst) worst = x;
		else if (x > best) best = x;
		delta = x - mean;
		mean = mean + (delta / (float)counter);
		var = var + (delta * (x - mean));
	}
	vector<float> output (5);
	output[0] = best;
	output[1] = worst;
	output[2] = mean; //the accuracy of mean is probably sufficient
	output[3] = var / (float)counter; // counter = POP_SIZE * NUM_ISLANDS
	output[4] = numZeros;
	return output;
}

/**
 * Computes the fitness of each individual within each island across the given test set.
 * Returns a two-dimensional vector of [islands * individuals][test instances]
 */
vector< vector<float> > calculatePhenotypeFitnessDetail(vector<TestInstance2>* testInstances, char classification, int startIsland, int endIsland) {
	if ((startIsland < 0) || (endIsland > NUM_ISLANDS)) { cout << "Error in calculateFitnessDetail; island index out of range " << startIsland << ", " << endIsland << endl; cout.flush(); exit(-1); }
	Individual2* tempIndividual;
	vector< vector<float> > diversetableH ((endIsland-startIsland)*POP_SIZE, vector<float>((*testInstances).size()));
	for (int i = startIsland; i < endIsland; ++i) {
		for (int j = 0; j < POP_SIZE; ++j) {
			tempIndividual = islands[i].getIndividual(j);
			(*tempIndividual).resetConfMat();
			for (int k = 0; k < (*testInstances).size(); ++k) {
				switch (classification) {
				case 'h': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = (*tempIndividual).classiHiFi(&(*testInstances)[k]); break;
				case 'l': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = (float)(*tempIndividual).classify(&(*testInstances)[k]); break;
				}
			}
		}
		// Not needed because the classify methods above do not change the individuals' fitness that is stored,
		// only the confmat that gets reset at the next finsess update in doOneGeneration
		//islands[i].a_pop->updatePopulationFitness(WHICH_CLASSIFY); //reset the fitness levels of the current island
	}
	return diversetableH;
}

/**
 * Add up the fitness values from that all the individuals for each test instance.
 * Returns a vector with as many items as numTestInstances, or the number of test instances.
 */
vector<float> calculatePhenotypeFitnessTotals(vector< vector<float> >* details, int numTestInstances) {
	vector<float> fitnessTotals (numTestInstances);
	for (int i = 0; i < numTestInstances; ++i) {
		fitnessTotals[i] = 0;
	}
	for (int i = 0; i < (*details).size(); ++i) {
		for (int j = 0; j < numTestInstances; ++j) {
			fitnessTotals[j] += (*details)[i][j];
		}
	}
	return fitnessTotals;
}

/**
 * Ranks values in the given vector starting with 1 as the largest value up to n as the smallest value.
 * Ranking values are not repeated.
 * Returns a vector of int that is the same size as the input vector.
 */
vector<int> calculateRankings(vector<float> fitnessTotals, bool usePriorityRanking) {
	vector<int> rankings (fitnessTotals.size());
	//int numerator = 0;
	for (int i = 0; i < fitnessTotals.size(); ++i) { // for every value, one at a time
		rankings[i] = 0;
		for (int j = 0; j < fitnessTotals.size(); ++j) {
			if ((fitnessTotals[i] <= fitnessTotals[j]) && ((i <= j) || (fitnessTotals[i] != fitnessTotals[j]))) { // i >= j makes sure we don't get duplicates in the rankings
				rankings[i] = rankings[i]+1;
			}
		}
		//if (fitnessTotals[i] == 0.0) numerator++;  // if nobody got this problem right, count it as non-diverse
	}
	if (usePriorityRanking) { // convects the ranking into a priority rank so that the vector is a set of indexes into a vector of weights according to their value.
		vector<int> priorityRank (fitnessTotals.size());
		for (int i = 0; i < fitnessTotals.size(); ++i) {
			priorityRank[(rankings[i]-1)] = i;
		}
		rankings = priorityRank;
	}
	return rankings;
}

/**
 * Scales the values such that their combined total equals the given target total.
 */
void scaleToTotal(vector<float>* vec, float targetTotal) {
	float scale = 0.0;
	//vector<float> output (vec.size());
	for (int i = 0; i < (*vec).size(); ++i) {
		scale += (*vec)[i];
	}
	scale = targetTotal / scale;
	for (int i = 0; i < (*vec).size(); ++i) {
		(*vec)[i] *= scale;
	}
}

/**
 * Scale the values in the given vector so that their max is the given limit and the other values are scaled accordingly.
 */
void scaleToLimit(vector<float>* vec, float limit) {
	float max = 0.0;
	//vector<float> output (vec.size());
	for (int i = 0; i < (*vec).size(); ++i) {
		if (max < (*vec)[i]) max = (*vec)[i];
	}
	float scale = limit / max;
	for (int i = 0; i < (*vec).size(); ++i) {
		(*vec)[i] *= scale;
	}
}

void scaleProportionally(vector<float>* vec, float scale) {
	for (int i = 0; i < (*vec).size(); ++i) {
		(*vec)[i] /= scale;
	}
}

/**
 * Returns a vector of weights where the largest number had the smallest total fitness in the given input vector of fitness totals.
 * Larger fitness values carry less weight and have smaller numbers in the output vector.
 * Returned vector is the same size as the input vector.
 */
vector<float> calculatePhenotypeFitnessWeights(vector<float>* fitnessTotals) {
	vector<float> weights ((*fitnessTotals).size());
	/*
	float total = 0.0;
	for (int i = 0; i < (*fitnessTotals).size(); ++i) { // for every value, one at a time
		weights[i] = 0;
		total += (*fitnessTotals)[i];
	}
	total = sqrt(total);*/
	for (int i = 0; i < (*fitnessTotals).size(); ++i) {
		weights[i] = ((float)POP_SIZE - (*fitnessTotals)[i]) / (float)POP_SIZE;
		//if ((*fitnessTotals)[i] != 0.0) weights[i] = total / (*fitnessTotals)[i];
		//else weights[i] = 0.0;  // what value should we ascribe to a test instance when nobody gets it right?
	}
	return weights;
}

/**
 * Returns a vector that represents the relevance of each individual.
 * Each individual's fitness on each test instance is multiplied by the weight associated with that test instance.
 * Returned vector is the same size as the input detaiedFitness vector or (islands * population-size).
 */
vector<float> calculateIndividualPhenotypeRelevance(vector< vector<float> >* detailedFitness, vector<float>* weights) {
	vector<float> relevanceVec ((*detailedFitness).size());
	if ((*detailedFitness)[0].size() != (*weights).size()) { // check the inputs
		cout << "Error at calculateIndividualRelevance: indiviual fitness vectors and fitness weights vector do not match" << endl;
		exit(-1);
	}
	for (int i = 0; i < (*detailedFitness).size(); ++i) {
		float individualValue = 0.0;
		for (int j = 0; j < (*weights).size(); ++j) {
			individualValue += (*detailedFitness)[i][j] * (*weights)[j]; // Sum the product of the individual's performance and the weight of each test instance.
		}
		relevanceVec[i] = individualValue;
	}
	return relevanceVec;
}

void getRelevanceByIsland(vector<float>* detail, vector <vector<float> >* islandRelevance) {
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		for (int j = 0; j < POP_SIZE; ++j) {
			(*islandRelevance)[i][j] = (*detail)[(POP_SIZE * i) + j];
		}
	}
}

/**
 * This wrapper function constructs the phenotype diversity of the population.
 */
vector<float> getPhenotypeRelevance(vector<TestInstance2>* testInstances, int startIsland, int endIsland, bool isOverallDiversity, vector <vector<float> >* islandRelevance, float relevanceWeight) {
	vector< vector<float> > detailedFitness = calculatePhenotypeFitnessDetail(testInstances, WHICH_CLASSIFY, startIsland, endIsland);  // this stores data in the array and returns detailed results also
	vector<float> fitnessTotals = calculatePhenotypeFitnessTotals(&detailedFitness, (*testInstances).size());
	vector<float> weights =  calculatePhenotypeFitnessWeights(&fitnessTotals);
	int numZeros = 0;
	for (int i = 0; i < weights.size(); ++i) {
		if (weights[i] == 1.0) numZeros++;  // checks to see how many test instances have 0 weight.
		if ((weights[i] < 0) || (weights[i] > 1.0)) {cout << endl << "error in phenotype weights: invalid value: " << weights[i] << endl; exit(-1);}
	}
	vector<float> individualRelevance = calculateIndividualPhenotypeRelevance(&detailedFitness, &weights);
	vector<float> pDiversity = getDiversityValues(&individualRelevance, numZeros);
	if (isOverallDiversity && (WHICH_SELECT == 3)) {
		// individualRelevance has values that have been accumulated over the number of test instances.
		// Therefore we use the number of test instances as the scale.
		scaleProportionally(&individualRelevance, (*testInstances).size());
		/*float max = 0.0;
		for (int i = 0; i < individualRelevance.size(); ++i) { if (individualRelevance[i] > max) max = individualRelevance[i]; }
		cout << "m: " << max << endl;*/
		if (relevanceWeight < 1.0) combineRelevanceWithFitness(relevanceWeight, &individualRelevance);
		getRelevanceByIsland(&individualRelevance, islandRelevance);
	}
	return pDiversity;
}

vector<float> getHammingRelevance(vector <vector<float> >* islandRelevance, bool isOverallDiversity, float relevanceWeight) {
	vector<float> relevance = getHammingRelevanceDetail();
	vector<float> diversityValues = getDiversityValues(&relevance, 0);
	if (isOverallDiversity && (WHICH_SELECT == 4)) {
		// relevance has values that have been accumulated over RULE_LEN.
		// Therefore we scale by the rule length.
		scaleProportionally(&relevance, RULE_LEN);
		/*float max = 0.0;
		for (int i = 0; i < relevance.size(); ++i) { if (relevance[i] > max) max = relevance[i]; }
		cout << "m: " << max << endl;*/
		if (relevanceWeight < 1.0) combineRelevanceWithFitness(relevanceWeight, &relevance);
		getRelevanceByIsland(&relevance, islandRelevance);
	}
	return diversityValues;
}

vector<TestInstance2> getTIVector(TestSet set, int num) {
	vector<TestInstance2> output (num);
	for (int i = 0; i < num; ++i) {
	  output[i] = * set.getTI(i);
	}
	return output;
}

void usage(){
	printf("Arguments to Singlenode version of GAchess executable:\n");
	printf(" 1.  test instance data file\n");
	printf(" 2.  logfile tag\n");
	printf(" 3.  migration strategy (r/s/w/f/n)\n");
	printf(" 4.  migration interval\n"); 
	printf(" 5.  num migrants per island\n");
	printf(" 6.  number of islands\n");
	printf(" 7.  num neighbors (must be less than the number of islands)\n");
	printf(" 8.  size of population (on each island)\n");
	printf(" 9.  num test cases to sample\n");
	printf(" 10. number of generations to simulate\n");
	printf(" 11. which fitness function to use (h/l)\n");
	printf(" 12. which classification function to use on full testset (h/l)\n");
	printf(" 13. (optional) mutation probability that given individual will have single-point mutation\n");
	printf(" 14. (optional) initial seed to select sample of test instances\n");
	printf(" 15. (optional) which class of test instances to use (only with 18 or fewer islands)\n");
}

int main(int argc, char * argv[])
/* args:
   argv[1] datafile,
   argv[2] logfile-tag,
   argv[3] migration-strategy (r:Rand, s:Strong, w:Weak, f:Fake, n:None),
   argv[4] migration interval (WHEN_MIGRATE),
   argv[5] number migrants per island,
   argv[6] number of islands,
   argv[7] number of neighbors,
   argv[8] size of population,
   argv[9] number of test cases to use,
   argv[10] number of generations to simulate,
   argv[11] which fitness function to use,
   argv[12] which classification function to use,
   argv[13] (optional) mutation probability that given individual will have single-point mutation
   argv[14] (optional) initial seed to use to select sample of test instances
   argv[15] (optional) depth of testinstance selecion (-1:draw, 0:zero, 1:one,...,16:sixteen)
 */
{
  srand( time(NULL) );
  int WHEN_PRINT_DATA = 1;
  int INIT_SEED;
  vector<TestInstance2> all_tests;

  //printf("SingleNode: Ready to start.\n");

  if (argc < 13) { usage(); exit(-1); }

  time_t start, end;

  // read filename from which to get test instances
  string filename = argv[1];

  WHICH_MIGRATION = argv[3][0];
  WHEN_MIGRATE = atoi(argv[4]);//atoi(argv[4]);
  NUM_MIGRANTS_PER_ISLAND = atoi(argv[5]);
  NUM_ISLANDS = atoi(argv[6]);
  NUM_NEIGHBORS = atoi(argv[7]);
  if (NUM_NEIGHBORS >= NUM_ISLANDS) { printf("Max number of neighbors is one less than number of islands\n"); usage(); exit(-1); }
  if ((WHICH_MIGRATION == 'n') && (NUM_NEIGHBORS != 0)) { usage(); exit(-1); }
  POP_SIZE = atoi(argv[8]);
  if (PROB_REMAIN * POP_SIZE < NUM_NEIGHBORS * NUM_MIGRANTS_PER_ISLAND) {printf("too many migrants per neighbor for given population size\n"); usage(); exit(-1); }
  NUM_IMMIGRANTS = NUM_MIGRANTS_PER_ISLAND * NUM_NEIGHBORS;
  NUM_TEST_CASES_TO_USE = atoi(argv[9]);
  MAX_GENERATIONS = atoi(argv[10]);
  WHICH_FITNESS = argv[11][0];
  WHICH_CLASSIFY = argv[12][0];
  if (argc >= 14)
    PROB_MUTATE = (int)(RAND_MAX * atof(argv[13]));
  else
    PROB_MUTATE = (int)(RAND_MAX * 1.0/POP_SIZE); // default: one individual in population has one bit flipped
  if (argc >= 15)
    INIT_SEED = atoi(argv[14]);
  else
    INIT_SEED = time(NULL) + 93108;
  //printf("Seed value: %d\n", INIT_SEED);
  //  DEPTH_OF_TEST_INSTANCES = atoi(argv[11]);

  //printf("SingleNode: about to open test instance file\n");
  //get test instances from file
  ifstream file;
  //cout << "l196" << endl;
  file.open(filename.c_str());
  string line;
  while(!file.eof())
    {
      //cout << "l201" << endl;
      getline(file, line);
      if (line != "")
	{
	  TestInstance2 ti = TestInstance2(line);
	  all_tests.push_back(ti);
	}
    }
  //cout << "l209" << endl;
  if (argc == 15 && NUM_ISLANDS > 18)
    { printf("Invalid combination of depth and number of islands\n"); usage(); exit(-1); }
  float currentWeight = (WHICH_SELECT < 3) ? RELEVANCE_END : RELEVANCE_START;
  while (currentWeight >= RELEVANCE_END ) {
    for (int cycle = 0; cycle < NUM_CYCLES; cycle++) {
      start = time(NULL);
      islands = new SingleNode[NUM_ISLANDS];
      //initialize all islands
      //printf("Initializing the islands (and populations, etc.)\n");
      TestSet ts;
      //cout << "l214" << endl;
      if (argc > 15)
	{
	  printf("Not yet supporting stratified islands with split training/test sets\n"); exit(-1);
	  for(int j=0; j<NUM_ISLANDS; j++)
	    {
	      ts = TestSet(&all_tests, NUM_TEST_CASES_TO_USE, j-1);
	      islands[j] = SingleNode(j, ts);
	      islands[j].updateNodeIntRules();
	    }
	}
      else
	{
	  // printf("Creating the TestSet from the total testinstances in the file\n");
	  ts = TestSet(&all_tests, NUM_TEST_CASES_TO_USE, (float)NUM_TEST_CASES_TO_USE/all_tests.size());
	  //printf("Finished creating the Testset\n");
	  for(int j=0; j<NUM_ISLANDS; j++) {
	    islands[j] = SingleNode(j, ts);
	    islands[j].updateNodeIntRules();
	  }
	}
      vector<TestInstance2> tsVector = getTIVector(ts, NUM_TEST_CASES_TO_USE);
      //initialize the logFile
      string s, outFile1, outFile2;
      stringstream out1, out2;

      bool depthdivision;
      if (argc == 15)
	{depthdivision = true;}
      else
	{depthdivision = false;}

      outFile1 += "log/sim.overall.";
      outFile2 += "log/sim.islands.";

      out1 << "rep-" << argv[2] << WHICH_FITNESS << "fifit-" << WHICH_CLASSIFY << "ficlsfy-" \
	   << WHICH_MIGRATION << "ms" \
	   << WHEN_MIGRATE << "mi" \
	   << NUM_MIGRANTS_PER_ISLAND << "mn" \
	   << NUM_ISLANDS << "i" \
	   << NUM_NEIGHBORS << "n" \
	   << POP_SIZE << "p" \
	   << NUM_TEST_CASES_TO_USE << "ti" \
	   << MAX_GENERATIONS << "g";
      out2 << "selection-" << WHICH_SELECT << "_relevance-" << \
	currentWeight << "_cycle-" << \
	cycle << ".prn";
      s = out2.str();
      outFile1 += s;
      outFile2 += s;
      s = out1.str();
      remove(outFile1.c_str());
      logFile1 = fopen(outFile1.c_str(), "w");
      remove(outFile2.c_str());
      logFile2 = fopen(outFile2.c_str(), "w");

      //run the generations
      //printf("Starting run with %d neighbors, %d migrants per island, %c migration strategy\n", NUM_NEIGHBORS, NUM_MIGRANTS_PER_ISLAND, WHICH_MIGRATION);
      //SingleNode i;
      float mostFit=0;
      int mostFitIsland=-1;

      // Alternate log file arrangement:
      fprintf(logFile1, "Log file for Overall Diversity\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tEVALUATED OVER ALL TESTS\nMost Fit\ton island\tat generation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\tTI with no weight\tMax Hamming-Estimated Diversity\tMin Hamming-Estimated Diversity\tAverage Hamming-Estimated Diversity\tHamming-Estimated Diversity Variance\t\tMost Fit\ton island\tat generation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\tTI with no weight");
      fprintf(logFile2, "Log file for Island Diversity\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tFull Testset\nMost Fit\tAverage Fitness\tof generation\ton island\tStandard Deviation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\tTI with no weight\t\tMost Fit\tAverage Fitness\tof generation\ton island\tStandard Deviation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\tTI with no weight\n");
      // allocating memory
      switch (WHICH_SELECT)
	{	case 0: cout << "Using random selection with no bias toward fitness or diversity." << endl; break;
	case 1: cout << "Using tournament fitness selection." << endl; break;
	case 2: cout << "Using alt fitness selection." << endl; break;
	case 3: cout << "Using phenotype relevance selection.\nThe selection bias toward relevance is " << currentWeight << "." << endl; break;
	case 4: cout << "Using hamming-estimated relevance selection.\nThe selection bias toward relevance is " << currentWeight << "." << endl; break;
	}
      vector <vector<float> > islandRelevance (NUM_ISLANDS, vector<float> (POP_SIZE));
      vector<float> HamDiv (4, 0.0);//getPairwiseHammingDiversityForAllNodes();
      vector<float> FitDiv = getFitnessDiversity(0, NUM_ISLANDS);
      vector<float> PhenDiv = getPhenotypeRelevance(&tsVector, 0, NUM_ISLANDS, true, &islandRelevance, currentWeight);
      vector<float> HamRel = getHammingRelevance(&islandRelevance, true, currentWeight);
      for (int node = 0; node < NUM_ISLANDS; ++node) {
	islands[node].updateNodeRelevance(&islandRelevance[node]);  // update the relevance at the first generation
      }
      fprintf(logFile1, "\nSTARTING POINT\t\t\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%f\t%f\t%f\t%f",
	      (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4], HamRel[0], HamRel[1], HamRel[2], HamRel[3]);
      cout << "Starting cycle " << cycle+1 << ". Using relevance weight of " << currentWeight << ".\nCounting generations: ";
      for(int gen=0; gen < MAX_GENERATIONS; gen++)
	{
	  cout << gen << " ";
	  cout.flush();
	  mostFit = 0.0;
	  for(int island=0; island < NUM_ISLANDS; island++) {
	    //i = islands[island];
	    islands[island].doOneGeneration(gen);
	    if ( (gen+1) % WHEN_PRINT_DATA == 0 ) {
	      if (islands[island].a_pop->getPopulationMaxFitness() > mostFit) {
		mostFit = islands[island].a_pop->getPopulationMaxFitness();
		mostFitIsland = island;
	      }
	      /*
		HamDiv = i.a_pop->getInternalHammingDiversity();
		FitDiv = getFitnessDiversity(island, (island+1));
		PhenDiv = getPhenotypeRelevance(tsVector, island, (island+1), false, &islandRelevance, currentWeight);
		fprintf(logFile2, "%f\t%f\t%i\t%i\t%f\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\n",
		i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4]);
	      */
	      fprintf(logFile2, "Most Fit: %f Average Fitness: %f of generation %i on island %i Standard Deviation: %f\n",
		      islands[island].a_pop->getPopulationMaxFitness(), islands[island].a_pop->getPopulationAvgFitness(), islands[island].a_pop->getGeneration(), island, islands[island].a_pop->getStdev());
	    }
	  }
	  //HamDiv = getPairwiseHammingDiversityForAllNodes();
	  FitDiv = getFitnessDiversity(0, NUM_ISLANDS);
	  PhenDiv = getPhenotypeRelevance(&tsVector, 0, NUM_ISLANDS, true, &islandRelevance, currentWeight);  // we use the training set untel it's time to use the full test set
	  HamRel = getHammingRelevance(&islandRelevance, true, currentWeight);
	  fprintf(logFile1, "\n%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%f\t%f\t%f\t%f",
		  mostFit, mostFitIsland, islands[0].a_pop->getGeneration(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4], HamRel[0], HamRel[1], HamRel[2], HamRel[3]);
	  for (int node = 0; node < NUM_ISLANDS; ++node) {
	    islands[node].updateNodeRelevance(&islandRelevance[node]);
	  }
	  if (((gen+1) % WHEN_FULL_TEST) == 0) {// need the +1 because population's gen-count has been updated during doOneGen
	    mostFit = 0;
	    for(int island=0; island < NUM_ISLANDS; island++)
	      {
		//i = islands[island];
		islands[island].a_pop->updatePopulationFitness(&all_tests,WHICH_CLASSIFY);
		if(islands[island].a_pop->getPopulationMaxFitness() > mostFit)
		  {
		    mostFit = islands[island].a_pop->getPopulationMaxFitness();
		    mostFitIsland = island;
		  }
		fprintf(logFile2, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%f\t%f\t%i\t%i\t%f\n",
			islands[island].a_pop->getPopulationMaxFitness(), islands[island].a_pop->getPopulationAvgFitness(), islands[island].a_pop->getGeneration(), island, islands[island].a_pop->getStdev());
		//i.a_pop->getBestIndividual().dumpConfMat(logFile2);
	      }
	    PhenDiv = getPhenotypeRelevance(&all_tests, 0, NUM_ISLANDS, false, &islandRelevance, currentWeight);
	    fprintf(logFile1, "\t\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i",
		    mostFit, mostFitIsland, islands[0].a_pop->getGeneration(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4]);
	    for (int x = 0; x < NUM_ISLANDS; ++x) {
	      islands[x].a_pop->updatePopulationFitness(WHICH_CLASSIFY); //reset the fitness levels in the islands before returning to training set
	    }
	  }
	  //if ( gen+1 > WHEN_FULL_TEST * 10 ) WHEN_FULL_TEST *= 10;
	}
      delete [] islands;
      cout << "\n\n";
      // print end time
      end = time(NULL);
      fprintf(logFile1, "\n\nElapsed time %f seconds or %f minutes\n%s", difftime(end, start), difftime(end, start)/60.0, s.c_str() );
      fprintf(logFile2, "\n\nElapsed time %f seconds or %f minutes\n%s", difftime(end, start), difftime(end, start)/60.0, s.c_str() );
      fclose(logFile1);
      fclose(logFile2);
    }
    currentWeight -= RELEVANCE_INCREMENT;
  }
  return 0;
}

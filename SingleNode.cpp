//#include "krktestinstance.h"
//#include "dataset.h"
#include "population.h"
#include "SingleNode.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

// GLOBAL VARIABLES
FILE *logFile; // for log data

char WHICH_MIGRATION; // r: Random, s: Strong, w: Weak, n: None, f: Fake random migrant, . . .
int WHEN_MIGRATE;
int NUM_MIGRANTS_PER_ISLAND;
int NUM_IMMIGRANTS;
int NUM_ISLANDS;
int NUM_NEIGHBORS;
int POP_SIZE;
int NUM_TEST_CASES_TO_USE;
int DEPTH_OF_TEST_INSTANCES;
int MAX_GENERATIONS;
char WHICH_FITNESS;
char WHICH_CLASSIFY;
int PROB_MUTATE; // an int < RAND_MAX representing probability a single individual will get "hit"
// use RULE_LEN instead: int genomeLength = NUM_FEATURES*RULE_CASES*8;

SingleNode * islands;

//template <class T>
SingleNode::SingleNode(int r, DataSet<KRKTestInstance> ts)
{
	myrank = r;
	//  printf("SingleNode: construct an instance of a single-node\n");
	//set the default values
	//  printf("SingleNode: about to create a population\n");
	a_pop = new Population(r, NUM_ISLANDS, ts, POP_SIZE);
	a_pop->updatePopulationIntRules();
	//  printf("SingleNode: have population and now update that pop's fitness\n");
	a_pop->updatePopulationFitness(WHICH_FITNESS, 't');
	customs = new Individual2[NUM_IMMIGRANTS];
}

//template <class T>
SingleNode::SingleNode()
// Default constructor with empty test-set
{

}

//template <class T>
void SingleNode::doOneGeneration(int thisgen, char whichSet)
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
	a_pop->updatePopulationFitness(WHICH_FITNESS, whichSet);
}

//template <class T>
void SingleNode::updateNodeIntRules() {
	a_pop->updatePopulationIntRules();
}

//template <class T>
void SingleNode::updateNodeRelevance(vector<float>* islandRelevance) {
	a_pop->updatePopulationRelevance(islandRelevance);
}

//template <class T>
Individual2* SingleNode::getIndividual(int index) {
	return a_pop->getIndividual(index);
}

//template <class T>
void SingleNode::addIslandBitTotal(float * totals) {
	a_pop->addPopulationBitTotal(totals);
}

//template <class T>
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
vector< vector<float> > calculatePhenotypeFitnessDetail(DataSet<KRKTestInstance>* testInstances, char classification, char whichSet, int startIsland, int endIsland) {
	if ((startIsland < 0) || (endIsland > NUM_ISLANDS)) { cout << "Error in calculateFitnessDetail; island index out of range " << startIsland << ", " << endIsland << endl; cout.flush(); exit(-1); }
	int numTestInstances = testInstances->trainSetSize();
	if (whichSet == 'f') numTestInstances = testInstances->testSetSize() ;
	Individual2* tempIndividual;
	KRKTestInstance tempTI;
	vector< vector<float> > diversetableH ((endIsland-startIsland)*POP_SIZE, vector<float>(numTestInstances));
	for (int i = startIsland; i < endIsland; ++i) {
		for (int j = 0; j < POP_SIZE; ++j) {
			tempIndividual = islands[i].getIndividual(j);
			(*tempIndividual).resetConfMat();
			for (int k = 0; k < numTestInstances; ++k) {
				if (whichSet == 't') tempTI = (KRKTestInstance)(*testInstances->getTrainI(k));
				else if (whichSet == 'f') tempTI = (KRKTestInstance)(*testInstances->getTestI(k));
				switch (classification) {
				case 'h': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = tempTI.fitnessHiFi(tempIndividual); break;
				case 'l': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = (float)(tempTI.classify(tempIndividual)); break;
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
		//if (((*fitnessTotals)[i] < 0) || ((*fitnessTotals)[i] > POP_SIZE)) {cout << "\ncalculatePhenotypeFitnessWeights: error in computing fitness weights. Value is " << (*fitnessTotals)[i] << ". Max is " << POP_SIZE << ".\n"; cout.flush(); exit(-1);}
		if ((*fitnessTotals)[i] > POP_SIZE) {
			//TODO : This is just a stop gap measure!
			(*fitnessTotals)[i] = POP_SIZE;
		}
		//TODO : Is fitnessTotals allowing there to be nubmers greater than 50? Is it a float precision issue?
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
		cout << "Error at calculateIndividualRelevance: indiviual fitness vectors (" << (*detailedFitness)[0].size() << ") and fitness weights vector (" << (*weights).size() << ") do not match" << endl;
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
vector<float> getPhenotypeRelevance(DataSet<KRKTestInstance>* testInstances, int startIsland, int endIsland, bool isOverallDiversity, bool isFullTestSet, vector <vector<float> >* islandRelevance, float relevanceWeight, char classification) {
	int numTestInstances = testInstances->trainSetSize();
	char whichSet = 't';
	if (isFullTestSet) {
		whichSet = 'f';
		numTestInstances = testInstances->testSetSize();
	}
	vector< vector<float> > detailedFitness = calculatePhenotypeFitnessDetail(testInstances, classification, whichSet, startIsland, endIsland);  // this stores data in the array and returns detailed results also
	vector<float> fitnessTotals = calculatePhenotypeFitnessTotals(&detailedFitness, numTestInstances);
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
		scaleProportionally(&individualRelevance, numTestInstances);
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

vector<KRKTestInstance*> getTIVector(DataSet<KRKTestInstance> * set, char whichSet) {
	int num = set->testSetSize();
	if (whichSet == 't') num = set->trainSetSize();
	vector<KRKTestInstance*> output (num);
	for (int i = 0; i < num; ++i) {
	  if (whichSet == 'f') output[i] = (KRKTestInstance*)set->getTestI(i);
	  else if (whichSet == 't') output[i] = (KRKTestInstance*)set->getTrainI(i);
	}
	return output;
}

void usage(){
	printf("Arguments to Singlenode version of GAchess executable:\n");
	printf(" 1.  test instance data file\n");
	printf(" 2.  logfile tag\n");
	printf(" 3.  migration strategy (r/s/w/f/n)\n");
	printf(" 4.  migration interval\n");
	printf(" 4. (optional) initial seed to select sample of test instances\n");
	printf(" 5.  num migrants per island\n");
	printf(" 6.  number of islands\n");
	printf(" 7.  num neighbors (must be less than the number of islands)\n");
	printf(" 8.  size of population (on each island)\n");
	printf(" 9.  num test cases to sample\n");
	printf(" 10. number of generations to simulate\n");
	printf(" 11. which fitness function to use (h/l)\n");
	printf(" 12. which classification function to use on full dataset (h/l)\n");
	printf(" 13. (optional) mutation probability that given individual will have single-point mutation\n");
	//printf(" 15. (optional) which class of test instances to use (only with 18 or fewer islands)\n");
}

void printDS(DataSet<KRKTestInstance>* ds) {
	KRKTestInstance * testI;
	int depth = ds->getTestI(0)->getDepth();
	cout << "\n" << "Test Set:" << endl;
	for (int i = 0; i < ds->testSetSize(); i++) {
		testI = ds->getTestI(i);
		if (depth == testI->getDepth()) cout << i << ":" << testI->getDepth() << ",  ";
	}
	depth = ds->getTrainI(0)->getDepth();
	cout << "\n\n" << "Train Set:" << endl;
	for (int i = 0; i < ds->trainSetSize(); i++) {
		testI = ds->getTrainI(i);
		if (depth == testI->getDepth()) cout << i << ":" << testI->getDepth() << ",  ";
	}
}

void printTIVector(vector<KRKTestInstance>* vec) {
	int depth = (*vec)[0].getDepth();
	int counts[20] = {0};
	cout << "\n" << "Complete vector Set within SingleNode:" << endl;
	for (int i = 0; i < vec->size(); i++) {
		//if (depth != (*vec)[i]->getDepth()) cout << i << ":" << (*vec)[i]->getDepth() << ",  ";
		counts[(*vec)[i].getDepth()]++;
	}
	cout << "Vector: " << (*vec).size() << ", " << depth << endl;
	for (int i = 0; i < 20; i++) {
		cout << i << ":" << counts[i] << endl;
	}
}


int main(int argc, char * argv[])
/* args:
   argv[1] datafile,
   argv[2] logfile-tag,
   argv[3] migration-strategy (r:Rand, s:Strong, w:Weak, f:Fake, n:None),
   argv[3] (optional) initial seed to use to select sample of test instances
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
   	   not used argv[15] (optional) depth of krktestinstance selecion (-1:draw, 0:zero, 1:one,...,16:sixteen)
 */
{
  //argc[3] can be EITHER a migration strategy or a random seed if there are no further args.
  long seed = RANDOM_SEED;
  if (argc == 4) {
	  seed = atol(argv[3]);
  }
  if (seed == 0) seed = time(NULL);
  cout << "Using Random seed of " << seed << endl;
  srand( seed );
  char trainingSet = 't';
  char testSet = 'f';
  int WHEN_PRINT_DATA = 1;
  vector<KRKTestInstance> all_tests;
  time_t start, end;
  string filename = "insData";
  string logFileName;;
  //if (argc < 13) { usage(); exit(-1); }
  if (argc > 1) {
    filename = argv[1]; // data file
  }
  if (argc > 2) {
	logFileName = argv[2];// log file tag
  }
  else {
	  stringstream ss;
	  ss << "log-" << seed;
	  logFileName = ss.str();
  }
  if (argc > 4 && argc < 13) {cout << "Invalid parameters." << endl; exit(-1);}
  if (argc > 4) { // all other inputs
    WHICH_MIGRATION = argv[3][0];
    WHEN_MIGRATE = atoi(argv[4]);
    NUM_MIGRANTS_PER_ISLAND = atoi(argv[5]);
    NUM_ISLANDS = atoi(argv[6]);
    NUM_NEIGHBORS = atoi(argv[7]);
    POP_SIZE = atoi(argv[8]);
    NUM_TEST_CASES_TO_USE = atoi(argv[9]);
    MAX_GENERATIONS = atoi(argv[10]);
    WHICH_FITNESS = argv[11][0];
    WHICH_CLASSIFY = argv[12][0];
  }
  else {
	  cout << "Using default parameters." << endl;
	  WHICH_MIGRATION = WHICH_MIGRATION_DEFAULT;
	  WHEN_MIGRATE = WHEN_MIGRATE_DEFAULT;
	  NUM_MIGRANTS_PER_ISLAND = NUM_MIGRANTS_PER_ISLAND_DEFAULT;
	  NUM_ISLANDS = NUM_ISLANDS_DEFAULT;
	  NUM_NEIGHBORS = NUM_NEIGHBORS_DEFAULT;
	  POP_SIZE = POP_SIZE_DEFAULT;
	  NUM_TEST_CASES_TO_USE = NUM_TEST_CASES_TO_USE_DEFAULT;
	  MAX_GENERATIONS = MAX_GENERATIONS_DEFAULT;
	  WHICH_FITNESS = WHICH_FITNESS_DEFAULT;
	  WHICH_CLASSIFY = WHICH_CLASSIFY_DEFAULT;
    }
  NUM_IMMIGRANTS = NUM_MIGRANTS_PER_ISLAND * NUM_NEIGHBORS;
  if (NUM_NEIGHBORS >= NUM_ISLANDS) { printf("Max number of neighbors is one less than number of islands\n"); usage(); exit(-1); }
  if ((WHICH_MIGRATION == 'n') && (NUM_NEIGHBORS != 0)) { usage(); exit(-1); }
  if (PROB_REMAIN * POP_SIZE < NUM_NEIGHBORS * NUM_MIGRANTS_PER_ISLAND) {printf("too many migrants per neighbor for given population size\n"); usage(); exit(-1); }
  if (argc >= 14)
    PROB_MUTATE = (int)(RAND_MAX * atof(argv[13]));
  else
    PROB_MUTATE = (int)(RAND_MAX * 1.0/POP_SIZE); // default: one individual in population has one bit flipped
  //get test instances from file
  ifstream file;
  file.open(filename.c_str());
  string line;
  while(!file.eof())
    {
      //cout << "l201" << endl;
      getline(file, line);
      if (line != "")
	{
	  KRKTestInstance ti = KRKTestInstance(line);
	  all_tests.push_back(ti);
	}
    }
  //cout << "l209" << endl;
  if (argc == 15 && NUM_ISLANDS > 18)
    { printf("Invalid combination of depth and number of islands\n"); usage(); exit(-1); }
  //initialize the logFile
        string s, outFile;
        stringstream out1, out2;
        outFile += "log/sim.";
        out1 << "rep-" << logFileName << WHICH_FITNESS << "fifit-" << WHICH_CLASSIFY << "ficlsfy-" \
  	   << WHICH_MIGRATION << "ms" \
  	   << WHEN_MIGRATE << "mi" \
  	   << NUM_MIGRANTS_PER_ISLAND << "mn" \
  	   << NUM_ISLANDS << "i" \
  	   << NUM_NEIGHBORS << "n" \
  	   << POP_SIZE << "p" \
  	   << NUM_TEST_CASES_TO_USE << "ti" \
  	   << MAX_GENERATIONS << "g";
        out2 << "selection-" << WHICH_SELECT << "-rep-" << logFileName;
        s = out2.str();
        outFile += s;
        s = out1.str();
        remove(outFile.c_str());
        logFile = fopen(outFile.c_str(), "w");
  start = time(NULL);
  float currentWeight = (WHICH_SELECT < 3) ? RELEVANCE_END : RELEVANCE_START;
  while (currentWeight >= RELEVANCE_END ) {
    for (int cycle = 0; cycle < NUM_CYCLES; cycle++) {
      int WHEN_FULL_TEST = WHEN_FULL_DATA_START;
      int threshold = WHEN_FULL_DATA_THRESHOLD;
      islands = new SingleNode[NUM_ISLANDS];
      //initialize all islands
      //printf("Initializing the islands (and populations, etc.)\n");
      DataSet<KRKTestInstance> ts;
      //cout << "l214" << endl;
      if (argc > 15)
	{
	  printf("Not yet supporting stratified islands with split training/test sets\n"); exit(-1);
	  for(int j=0; j<NUM_ISLANDS; j++)
	    {
	      ts = DataSet<KRKTestInstance>(&all_tests, NUM_TEST_CASES_TO_USE, (float)j-1);
	      islands[j] = SingleNode(j, ts);
	      islands[j].updateNodeIntRules();
	    }
	}
      else
	{
	  // printf("Creating the DataSet from the total krktestinstances in the file\n");
	  ts = DataSet<KRKTestInstance>(&all_tests, NUM_TEST_CASES_TO_USE, (float)NUM_TEST_CASES_TO_USE/all_tests.size());
	  //printTIVector(&all_tests);
	  //printf("Finished creating the Testset\n");
	  for(int j=0; j<NUM_ISLANDS; j++) {
	    islands[j] = SingleNode(j, ts);
	    islands[j].updateNodeIntRules();
	  }
	}
      // we move this into DataSet and call based on the cirmstance.
      //vector<KRKTestInstance*> tsVector = getTIVector(&ts, NUM_TEST_CASES_TO_USE);

      bool depthdivision;
      if (argc == 15)
	{depthdivision = true;}
      else
	{depthdivision = false;}

    fprintf(logFile, "STARTING POINT\nStarting relevance %f on cycle %i\n", currentWeight, cycle);
      //run the generations
      //printf("Starting run with %d neighbors, %d migrants per island, %c migration strategy\n", NUM_NEIGHBORS, NUM_MIGRANTS_PER_ISLAND, WHICH_MIGRATION);
      float mostFit = 0.0;
      int mostFitIsland=-1;

    switch (WHICH_SELECT) {
      case 0: cout << "Using random selection with no bias toward fitness or diversity." << endl; break;
	  case 1: cout << "Using tournament fitness selection." << endl; break;
	  case 2: cout << "Using alt fitness selection." << endl; break;
	  case 3: cout << "Using phenotype relevance selection.\nThe selection bias toward relevance is " << currentWeight << "." << endl; break;
	  case 4: cout << "Using hamming-estimated relevance selection.\nThe selection bias toward relevance is " << currentWeight << "." << endl; break;
	}
      vector <vector<float> > islandRelevance (NUM_ISLANDS, vector<float> (POP_SIZE));
      vector<float> HamDiv (4, 0.0);//getPairwiseHammingDiversityForAllNodes();
      vector<float> FitDiv = getFitnessDiversity(0, NUM_ISLANDS);
      vector<float> PhenDiv = getPhenotypeRelevance(&ts, 0, NUM_ISLANDS, true, false, &islandRelevance, currentWeight, WHICH_FITNESS);
      vector<float> HamRel = getHammingRelevance(&islandRelevance, true, currentWeight);
      for (int node = 0; node < NUM_ISLANDS; ++node) {
	islands[node].updateNodeRelevance(&islandRelevance[node]);  // update the relevance at the first generation
      }
     fprintf(logFile, "Overall Most Fit %f on island %i at generation %i Max Hamming Diversity %i Min Hamming Diversity %i Average Hamming Diversity %f Hamming Diversity Variance %f Max Fitness Diversity %f Min Fitness Diversity %f Average Fitness Diversity %f Fitness Diversity Variance %f Max Phenotype Diversity %f Min Phenotype Diversity %f Average Phenotype Diversity %f Phenotype Diversity Variance %f TI with no weight %i Max Hamming-Estimated Diversity %f Min Hamming-Estimated Diversity %f Average Hamming-Estimated Diversity %f Hamming-Estimated Diversity Variance %f\n",
	      0.0, 0, 0, (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4], HamRel[0], HamRel[1], HamRel[2], HamRel[3]);
      cout << "Starting cycle " << cycle+1 << ". Using relevance weight of " << currentWeight << ".\nCounting generations: ";
      for(int gen=0; gen < MAX_GENERATIONS; gen++)
	{
	  cout << gen << " ";
	  cout.flush();
	  mostFit = 0.0;
	  for(int island=0; island < NUM_ISLANDS; island++) {
	    islands[island].doOneGeneration(gen,trainingSet);
	    if ( (gen+1) % WHEN_PRINT_DATA == 0 ) {
	      if (islands[island].a_pop->getPopulationMaxFitness() > mostFit) {
		mostFit = islands[island].a_pop->getPopulationMaxFitness();
		mostFitIsland = island;
	      }
	      /*
		HamDiv = i.a_pop->getInternalHammingDiversity();
		FitDiv = getFitnessDiversity(island, (island+1));
		PhenDiv = getPhenotypeRelevance(tsVector, island, (island+1), false, &islandRelevance, currentWeight);
	      */
	      fprintf(logFile, "Island Most Fit %f Average Fitness %f of generation %i on island %i Standard Deviation %f\n",
		      islands[island].a_pop->getPopulationMaxFitness(), islands[island].a_pop->getPopulationAvgFitness(), islands[island].a_pop->getGeneration(), island, islands[island].a_pop->getStdev());
	    }
	  }
	  //HamDiv = getPairwiseHammingDiversityForAllNodes(); // This operation in very expensive.
	  FitDiv = getFitnessDiversity(0, NUM_ISLANDS);
	  PhenDiv = getPhenotypeRelevance(&ts, 0, NUM_ISLANDS, true, false, &islandRelevance, currentWeight, WHICH_FITNESS);  // we use the training set untel it's time to use the full test set
	  HamRel = getHammingRelevance(&islandRelevance, true, currentWeight);
	  fprintf(logFile, "Overall Most Fit %f on island %i at generation %i Max Hamming Diversity %i Min Hamming Diversity %i Average Hamming Diversity %f Hamming Diversity Variance %f Max Fitness Diversity %f Min Fitness Diversity %f Average Fitness Diversity %f Fitness Diversity Variance %f Max Phenotype Diversity %f Min Phenotype Diversity %f Average Phenotype Diversity %f Phenotype Diversity Variance %f TI with no weight %i Max Hamming-Estimated Diversity %f Min Hamming-Estimated Diversity %f Average Hamming-Estimated Diversity %f Hamming-Estimated Diversity Variance %f\n",
			  mostFit, mostFitIsland, islands[0].a_pop->getGeneration(),
			  (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3],
			  FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3],
			  PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4],
			  HamRel[0], HamRel[1], HamRel[2], HamRel[3]);
	  for (int node = 0; node < NUM_ISLANDS; ++node) {
	    islands[node].updateNodeRelevance(&islandRelevance[node]);
	  }
	  if (((gen+1) % WHEN_FULL_TEST) == 0) {// need the +1 because population's gen-count has been updated during doOneGen
	    if ((gen+1) >= threshold) {
	    	threshold *= WHEN_FULL_DATA_MULTIPLY;
	    	WHEN_FULL_TEST *= WHEN_FULL_DATA_MULTIPLY;
	    }
		mostFit = 0;
	    for(int island=0; island < NUM_ISLANDS; island++)
	      {
		  islands[island].a_pop->updatePopulationFitness(WHICH_CLASSIFY, testSet);
		  if(islands[island].a_pop->getPopulationMaxFitness() > mostFit)
		  {
		    mostFit = islands[island].a_pop->getPopulationMaxFitness();
		    mostFitIsland = island;
		  }
		fprintf(logFile, "FULL TEST SET Island Most Fit %f Average Fitness %f of generation %i on island %i Standard Deviation %f\n",
			islands[island].a_pop->getPopulationMaxFitness(), islands[island].a_pop->getPopulationAvgFitness(), islands[island].a_pop->getGeneration(), island, islands[island].a_pop->getStdev());
		//islands[island].a_pop->getBestIndividual().dumpConfMat(logFile);
	      }
	    PhenDiv = getPhenotypeRelevance(&ts, 0, NUM_ISLANDS, true, true, &islandRelevance, currentWeight, WHICH_CLASSIFY);
	    fprintf(logFile, "FULL TEST SET Overall Most Fit %f on island %i at generation %i Max Hamming Diversity %i Min Hamming Diversity %i Average Hamming Diversity %f Hamming Diversity Variance %f Max Fitness Diversity %f Min Fitness Diversity %f Average Fitness Diversity %f Fitness Diversity Variance %f Max Phenotype Diversity %f Min Phenotype Diversity %f Average Phenotype Diversity %f Phenotype Diversity Variance %f TI with no weight %i Max Hamming-Estimated Diversity %f Min Hamming-Estimated Diversity %f Average Hamming-Estimated Diversity %f Hamming-Estimated Diversity Variance %f\n",
		    mostFit, mostFitIsland, islands[0].a_pop->getGeneration(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3], (int)PhenDiv[4], HamRel[0], HamRel[1], HamRel[2], HamRel[3]);
	    for (int x = 0; x < NUM_ISLANDS; ++x) {
	      islands[x].a_pop->updatePopulationFitness(WHICH_FITNESS, trainingSet); //reset the fitness levels in the islands before returning to training set
	    }
	  }
	}
      delete [] islands;
      fprintf(logFile, "\n\n");
      cout << "\n\n";
    }
    currentWeight -= RELEVANCE_INCREMENT;
  }
  // print end time
  end = time(NULL);
  fprintf(logFile, "\n\nElapsed time %f seconds or %f minutes\n%s", difftime(end, start), difftime(end, start)/60.0, s.c_str() );
  fprintf(logFile,"\nParameters used:\nRandom Seed %ld\nWHICH_MIGRATION %c\nWHEN_MIGRATE %i\nNUM_MIGRANTS_PER_ISLAND %i\nNUM_ISLANDS %i\nNUM_NEIGHBORS %i\nPOP_SIZE %i\nNUM_TEST_CASES_TO_USE %i\nMAX_GENERATIONS %i\nWHICH_FITNESS %c\nWHICH_CLASSIFY %c\n",seed,WHICH_MIGRATION,WHEN_MIGRATE,NUM_MIGRANTS_PER_ISLAND,NUM_ISLANDS,NUM_NEIGHBORS,POP_SIZE,NUM_TEST_CASES_TO_USE,MAX_GENERATIONS,WHICH_FITNESS,WHICH_CLASSIFY);
  fprintf(logFile, "\n\nConfig.h contents:\nNUM_FEATURES %i\nINST_LEN %i\nN_ULONG_IN_RULE %i\nRULE_CASES %i\nRULE_LEN %i\nONE_FREQ %i\nMAX_NUM_CROSS_BREED %i\nPROB_REMAIN %f\nTAG %i\nWHICH_SELECT %i\nTOURNAMENT_SIZE %i\nRELEVANCE_START %f\nRELEVANCE_END %f\nRELEVANCE_INCREMENT %f\nNUM_CYCLES %i\nRANDOM_SEED %i\nWHEN_FULL_DATA_MULTIPLY %i\nWHEN_FULL_DATA_THRESHOLD %i\nWHEN_FULL_DATA_START %i\nWHICH_MIGRATION_DEFAULT %c\nWHEN_MIGRATE_DEFAULT %i\nNUM_MIGRANTS_PER_ISLAND_DEFAULT %i\nNUM_ISLANDS_DEFAULT %i\nNUM_NEIGHBORS_DEFAULT %i\nPOP_SIZE_DEFAULT %i\nNUM_TEST_CASES_TO_USE_DEFAULT %i\nMAX_GENERATIONS_DEFAULT %i\nWHICH_FITNESS_DEFAULT %c\nWHICH_CLASSIFY_DEFAULT %c\n",NUM_FEATURES,INST_LEN,N_ULONG_IN_RULE,RULE_CASES,RULE_LEN,ONE_FREQ,MAX_NUM_CROSS_BREED,PROB_REMAIN,TAG,WHICH_SELECT,TOURNAMENT_SIZE,RELEVANCE_START,RELEVANCE_END,RELEVANCE_INCREMENT,NUM_CYCLES,RANDOM_SEED,WHEN_FULL_DATA_MULTIPLY,WHEN_FULL_DATA_THRESHOLD,WHEN_FULL_DATA_START,WHICH_MIGRATION_DEFAULT,WHEN_MIGRATE_DEFAULT,NUM_MIGRANTS_PER_ISLAND_DEFAULT,NUM_ISLANDS_DEFAULT,NUM_NEIGHBORS_DEFAULT,POP_SIZE_DEFAULT,NUM_TEST_CASES_TO_USE_DEFAULT,MAX_GENERATIONS_DEFAULT,WHICH_FITNESS_DEFAULT,WHICH_CLASSIFY_DEFAULT);


  fclose(logFile);
  return 0;
}

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
/*int*/float WHEN_MIGRATE;
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

SingleNode * islands;

SingleNode::SingleNode(int r, TestSet ts)
{
  myrank = r;
  //  printf("SingleNode: construct an instance of a single-node\n");
  //set the default values

  //  printf("SingleNode: about to create a population\n");
  a_pop = new Population(r, NUM_ISLANDS, ts, POP_SIZE);
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

  float randmigrate;
  srand( time(NULL) );
  randmigrate = ((float)rand())/RAND_MAX; 
  // send
  //printf("Node %d: sendMigrantStrings()\n", myrank);
  if (WHICH_MIGRATION != 'n' &&  (randmigrate > WHEN_MIGRATE)/*!(thisgen%WHEN_MIGRATE)*/ )
    sendMigrants();
  
  // survive
  //printf("Node %d: a_pop->selectToSurvive(%d)\n", myrank, ((int)(POP_SIZE*PROB_REMAIN)));
  a_pop->selectToSurvive((int)(POP_SIZE*PROB_REMAIN));

  // receive
  //printf("Node %d: receiveMigrantStrings()\n", myrank);
  // ignore when only one island, but need to receive migrants (unless sendMigrants takes care of it)
  if (WHICH_MIGRATION != 'n' && (randmigrate > WHEN_MIGRATE)/*!(thisgen%WHEN_MIGRATE)*/ )
    a_pop->processImmigrants(islands[myrank].customs, NUM_IMMIGRANTS);

  // breed
  //printf("Node %d: a_pop->generateOffspring(%d)\n", myrank,(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + NUM_IMMIGRANTS)));
  a_pop->generateOffspring(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + ((WHICH_MIGRATION != 'n' && (randmigrate > WHEN_MIGRATE)/*!(thisgen%WHEN_MIGRATE)*/ ) ? NUM_IMMIGRANTS : 0)));
  //a_pop->generateOffspring(POP_SIZE - (((int)(POP_SIZE*PROB_REMAIN)) + NUM_IMMIGRANTS));

  // switch to next gen
  //printf("Node %d: a_pop->nextGeneration()\n", myrank);
  a_pop->nextGeneration(PROB_MUTATE);

  // update population fitness
  //printf("Node %d: a_pop-updatePopulationFitness()\n", myrank);
  a_pop->updatePopulationFitness(WHICH_FITNESS);
}
  
Individual2 SingleNode::getIndividual(int index) {
	return a_pop->getIndividual(index);
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
vector<float> getHammingDiversityForAllNodes() {
  // get number of nodes NUM_ISLANDS
  // get number of indivuals in each node POP_SIZE
  // nested for loop
  int TOTAL_POP = POP_SIZE * NUM_ISLANDS;
  int bestDiversity = 0;
  int worstDiversity = 1000; // Start with a high number.
  vector<int> tempV;
  Individual2 currentIndividual;
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

void seeInside(vector<float> input, string leadStr) {
	bool bad = false;
	cout << endl << leadStr;
	for (int i = 0; i < input.size(); ++i) {
		if ((i % 10) == 0) cout << endl;
		if (isnanf(input[i]) != 0) bad = true;
		cout << input[i] << ", ";
	}
	cout << endl;
	cout.flush();
	if (bad) {cout << endl << endl << "Found nan!" << endl; cout.flush(); exit(-1);}
}

/*
 * Computes the variance of the fitness of the individuals in the entire population with the given test set.
 * Returned vector is <best, worst, average, variance>
 */
vector<float> getFitnessDiversity(vector<TestInstance2> testInstances, char classification, int startIsland, int endIsland) {
	if ((startIsland < 0) || (endIsland > NUM_ISLANDS)) { cout << "Error in getFitnessDiversity; island index out of range " << startIsland << ", " << endIsland << endl; exit(-1); }
	int counter = 0;
	float worst = 1000.0;
	float mean, x, delta, var, best = 0.0;
	for (int i = startIsland; i < endIsland; ++i) {
		islands[i].a_pop->updatePopulationFitness(testInstances, classification);
		for (int j = 0; j < POP_SIZE; ++j) {
			++counter;
			x = islands[i].getIndividual(j).getFitness();
			if (x < worst) worst = x;
			else if (x > best) best = x;
			//average += x;
			delta = x - mean;
			mean = mean + (delta / (float)counter);
			var = var + (delta * (x - mean));
		}
	}
	//average = average / (float) counter;
	//if (counter != (POP_SIZE*NUM_ISLANDS)) cout << "counter: " << counter << " total islands: " << POP_SIZE * NUM_ISLANDS << endl; // would indicate a problem in the counter
	//if (mean != average) cout << "mean: " << mean << " Average: " << average << endl; // Very small floating point error in calculating mean
	//return var;
	vector<float> output (4);
	output[0] = best;
	output[1] = worst;
	output[2] = mean; //the accuracy of mean is probably sufficient
	output[3] = var / (float)counter; // counter = POP_SIZE * NUM_ISLANDS
	return output;
}

/*
 * Computes the max, min, average, and variance of the numbers in the given vector.
 * Returned vector is <best, worst, average, variance>
 */
vector<float> getDiversityValues(vector<float> values) {
	int counter = 0;
	float worst = 1000.0;
	float mean, x, delta, var, best = 0.0;
	for (int i = 0; i < values.size(); ++i) {
		++counter;
		x = values[i];
		if (x < worst) worst = x;
		else if (x > best) best = x;
		delta = x - mean;
		mean = mean + (delta / (float)counter);
		var = var + (delta * (x - mean));
	}
	vector<float> output (4);
	output[0] = best;
	output[1] = worst;
	output[2] = mean; //the accuracy of mean is probably sufficient
	output[3] = var / (float)counter; // counter = POP_SIZE * NUM_ISLANDS
	return output;
}

/**
 * Computes the fitness of each individual within each island across the given test set.
 * Returns a two-dimensional vector of [islands * individuals][test instances]
 */
vector< vector<float> > calculateFitnessDetail(vector<TestInstance2> testInstances, char classification, int startIsland, int endIsland) {
	if ((startIsland < 0) || (endIsland > NUM_ISLANDS)) { cout << "Error in calculateFitnessDetail; island index out of range " << startIsland << ", " << endIsland << endl; cout.flush(); exit(-1); }
	SingleNode tempIsland;
	Individual2 tempIndividual;
	vector< vector<float> > diversetableH ((endIsland-startIsland)*POP_SIZE, vector<float>(testInstances.size())); // for hi-fi classification
	//vector< vector< vector<int> > > diversetableL (NUM_ISLANDS, vector< vector<int> >(POP_SIZE, vector<int>(testInstances.size())));  // for low-fi classification, but we need to return something else
	for (int i = startIsland; i < endIsland; ++i) {
		tempIsland = islands[i];
		for (int j = 0; j < POP_SIZE; ++j) {
			tempIndividual = tempIsland.getIndividual(j);
			tempIndividual.resetConfMat();
			for (int k = 0; k < testInstances.size(); ++k) {
				switch (classification) {
				case 'h': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = tempIndividual.classiHiFi(testInstances[k]); break;
				case 'l': diversetableH[((i-startIsland)*POP_SIZE)+j][k] = (float)tempIndividual.classify(testInstances[k]); break;
				}
			}
		}
	}
	return diversetableH;
}

/**
 * Add up the fitness values from that all the individuals for each test instance.
 * Returns a vector with as many items as numTestInstances, or the number of test instances.
 */
vector<float> calculateFitnessTotals(vector< vector<float> > details, int numTestInstances) {
	vector<float> fitnessTotals (numTestInstances);
	for (int i = 0; i < numTestInstances; ++i) {
		fitnessTotals[i] = 0;
	}
	for (int i = 0; i < details.size(); ++i) {
		for (int j = 0; j < numTestInstances; ++j) {
			fitnessTotals[j] += details[i][j];
		}
	}
	return fitnessTotals;
}

/**
 * Ranks values in the given vector starting with 1 as the largest value up to n as the smallest value.
 * Ranking values are not repeated.
 * Returns a vector of int that is the same size as the input vector.
 */
vector<int> calculateFitnessRankings(vector<float> fitnessTotals) {
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
	//return (float)1.0-(numerator/(float)testInstances.size());
	return rankings;
}

/**
 * Scales the weights such that their combined total equals the given target total.
 */
vector<float> scaleToTotal(vector<float> vec, float targetTotal) {
	float scale = 0.0;
	vector<float> output (vec.size());
	for (int i = 0; i < vec.size(); ++i) {
		scale += vec[i];
	}
	scale = targetTotal / scale;
	for (int i = 0; i < vec.size(); ++i) {
		output[i] = vec[i] * scale;
	}
	return output;
}

/**
 * Scale the weights in the given vector so that their max is the given limit and the other values are scaled accordingly.
 */
vector<float> scaleToLimit(vector<float> vec, float limit) {
	float max = 0.0;
	vector<float> output (vec.size());
	for (int i = 0; i < vec.size(); ++i) {
		if (max < vec[i]) max = vec[i];
	}
	float scale = limit / max;
	for (int i = 0; i < vec.size(); ++i) {
		output[i] = vec[i] * scale;
	}
	return output;
}

/**
 * Returns a vector of weights where the largest number had the smallest total fitness in the given input vector of fitness totals.
 * Larger fitness values carry less weight and have smaller numbers in the output vector.
 * Returned vector is the same size as the input vector.
 */
vector<float> calculateFitnessWeights(vector<float> fitnessTotals) {
	vector<float> weights (fitnessTotals.size());
	float total = 0.0;
	for (int i = 0; i < fitnessTotals.size(); ++i) { // for every value, one at a time
		weights[i] = 0;
		total += fitnessTotals[i];
	}
	for (int i = 0; i < fitnessTotals.size(); ++i) {
		float temp = total / fitnessTotals[i];
		/*if (isnanf(temp*0.0) != 0){// != 0) { //  Debug: infinite values are returned if the denominator is 0
			cout << endl << "Found inf number at i:" << i << " and " <<  total << " / " << fitnessTotals[i] << endl;
			seeInside(fitnessTotals, "Fitness Totals that ended the game:");
			exit(-1);
		}*/
		weights[i] = temp;
	}
	return weights;
}

/**
 * Returns a vector that represents the relavance of each individual.
 * Each individual's fitness on each test instance is multiplied by the weight associated with that test instance.
 * Returned vector is the same size as the input detaiedFitness vector or (islands * population-size).
 */
vector<float> calculateIndividualRelavance(vector< vector<float> > detailedFitness, vector<float> weights) {
	vector<float> relavanceVec (detailedFitness.size());
	if (detailedFitness[0].size() != weights.size()) { // check the inputs
		cout << "Error at calculateIndividualRelavance: indiviual fitness vectors and fitness weights vector do not match" << endl;
		exit(-1);
	}
	for (int i = 0; i < detailedFitness.size(); ++i) {
		float individualValue = 0.0;
		for (int j = 0; j < weights.size(); ++j) {
			individualValue += detailedFitness[i][j] * weights[j]; // Sum the product of the individual's performance and the veight of each test instance.
			/*if (isnanf(individualValue) != 0) {
				cout << endl << "nan i:" << i << ", j:" << j << "= "<< detailedFitness[i][j] << "*" << weights[j] << " and " << "vector: " << endl;
				seeInside(weights, "weights that caused error: ");
				exit(-1); }*/
		}
		relavanceVec[i] = individualValue;
	}
	return relavanceVec;
}

/**
 * This wrapper function constructs the phenotype diversity of the population.
 */
vector<float> getPhenotypeDiversity(vector<TestInstance2> testInstances, char classification, int startIsland, int endIsland) {
	vector< vector<float> > detailedFitness = calculateFitnessDetail(testInstances, classification, startIsland, endIsland);  // this stores data in the array and returns detailed results also
	vector<float> fitnessTotals = calculateFitnessTotals(detailedFitness, testInstances.size());
	//vector<int> rankings =  calculateFitnessRankings(fitnessTotals);
	vector<float> weights =  calculateFitnessWeights(fitnessTotals);
	//weights = scaleToTotal(weights, 1.0); // Scales the weights such that their total equals the scale
	//weights = scaleToLimit(weights, 1.0); // Scale the weights so that their max is the given limit and the other values are scaled accordingly.
	vector<float> individualRelavance = calculateIndividualRelavance(detailedFitness, weights);
	//seeInside(fitnessTotals, "FitnessTotals");
	//seeInside(weights, "weights");
	//seeInside(individualRelavance, "individualRelavance");
	return getDiversityValues(individualRelavance);
}

vector<TestInstance2> getTIVector(TestSet set, int num) {
	vector<TestInstance2> output (num);
	for (int i = 0; i < num; ++i) {
		output[i] = set.getTI(i);
	}
	return output;
}

void usage(){
  printf("Arguments to Singlenode version of GAchess executable:\n");
  printf(" 1.  test instance data file\n");
  printf(" 2.  logfile tag\n");
  printf(" 3.  migration strategy (r/s/w/f/n)\n");
  printf(" 4.  migration threshold\n"); //"interval\n"); 
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
  int WHEN_PRINT_DATA = 1;
  int INIT_SEED;
  vector<TestInstance2> all_tests;  // vector of random samples from test instances

  //printf("SingleNode: Ready to start.\n");

  if (argc < 13) { usage(); exit(-1); }

  time_t start, end;
  start = time(NULL);

  // read filename from which to get test instances
  string filename = argv[1];

  WHICH_MIGRATION = argv[3][0];
  WHEN_MIGRATE = atof(argv[4]);//atoi(argv[4]);
  NUM_MIGRANTS_PER_ISLAND = atoi(argv[5]);
  NUM_ISLANDS = atoi(argv[6]);
  // do not need to check limit on number of islands in SingleNode
  islands = new SingleNode[NUM_ISLANDS];
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

  //initialize all islands
  //printf("Initializing the islands (and populations, etc.)\n");
  TestSet ts;
  //cout << "l214" << endl;
  if (argc > 15) 
    {
      printf("Not yet supporting stratified islands with split training/test sets\n"); exit(-1);
      for(int j=0; j<NUM_ISLANDS; j++)
	{
	  ts = TestSet(all_tests, NUM_TEST_CASES_TO_USE, j-1);
	  islands[j] = SingleNode(j, ts);
	}
    }
  else
    {
      // printf("Creating the TestSet from the total testinstances in the file\n");
      ts = TestSet(all_tests, NUM_TEST_CASES_TO_USE, (float)NUM_TEST_CASES_TO_USE/all_tests.size());
      //printf("Finished creating the Testset\n");
      for(int j=0; j<NUM_ISLANDS; j++)
	{
	  islands[j] = SingleNode(j, ts);
	}
    }
  vector<TestInstance2> tsVector = getTIVector(ts, NUM_TEST_CASES_TO_USE);
  //initialize the logFile
  string s, outFile1, outFile2;
  stringstream out;

  bool depthdivision;
  if (argc == 15)
    {depthdivision = true;}
  else 
    {depthdivision = false;}
  
  outFile1 += "log/sim.overall.";
  outFile2 += "log/sim.islands.";

  out << "rep-" << argv[2] << WHICH_FITNESS << "fifit-" << WHICH_CLASSIFY << "ficlsfy-" \
      << WHICH_MIGRATION << "ms" \
      << WHEN_MIGRATE << "mi" \
      << NUM_MIGRANTS_PER_ISLAND << "mn" \
      << NUM_ISLANDS << "i" \
      << NUM_NEIGHBORS << "n" \
      << POP_SIZE << "p" \
      << NUM_TEST_CASES_TO_USE << "ti" \
      << MAX_GENERATIONS << "g" << ".prn";  // use a file extension that can be opened as a spreadsheet
  s = out.str();
  outFile1 += s;
  outFile2 += s;

  remove(outFile1.c_str());
  logFile1 = fopen(outFile1.c_str(), "w");
  remove(outFile2.c_str());
    logFile2 = fopen(outFile2.c_str(), "w");

  //run the generations
  //printf("Starting run with %d neighbors, %d migrants per island, %c migration strategy\n", NUM_NEIGHBORS, NUM_MIGRANTS_PER_ISLAND, WHICH_MIGRATION);
  SingleNode i;
  float mostFit=0;
  int mostFitIsland=-1;

  cout << "Counting generations: ";
  // Alternate log file arrangement:
  fprintf(logFile1, "Log file for Overall Diversity\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tEVALUATED OVER ALL TESTS\nMost Fit\ton island\tat generation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\t\tMost Fit\ton island\tat generation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance");
  fprintf(logFile2, "Log file for Island Diversity\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tFull Testset\nMost Fit\tAverage Fitness\tof generation\ton island\tStandard Deviation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\t\tMost Fit\tAverage Fitness\tof generation\ton island\tStandard Deviation\tMax Hamming Diversity\tMin Hamming Diversity\tAverage Hamming Diversity\tHamming Diversity Variance\tMax Fitness Diversity\tMin Fitness Diversity\tAverage Fitness Diversity\tFitness Diversity Variance\tMax Phenotype Diversity\tMin Phenotype Diversity\tAverage Phenotype Diversity\tPhenotype Diversity Variance\n");

  for(int gen=0; gen < MAX_GENERATIONS; gen++)
    {
      cout << gen << " ";
      cout.flush();
      //printf("Starting generation %d\n", gen);
      mostFit = 0.0;
      // later, we'll add a loop over the islands that are being simulated on this one node
      for(int island=0; island < NUM_ISLANDS; island++)
	{
	  //printf("Island %d: about to try one generation\n",island);
	  // ***** Instead of i, try directly accessing islands[island].method(...) below, to check for speed up??
	  i = islands[island];

	  //printf("...........OK, try this generation(%d)\n", gen);
	  i.doOneGeneration(gen);
	  //printf("Node %d: ran generation %d\n", island, gen);

	  if( (gen+1) % WHEN_PRINT_DATA == 0 )
	    {
	      if(i.a_pop->getPopulationMaxFitness() > mostFit)
		{
		  mostFit = i.a_pop->getPopulationMaxFitness();
		  mostFitIsland = island;
		}
	      vector<float> HamDiv = i.a_pop->getInternalHammingDiversity();
	      vector<float> FitDiv = getFitnessDiversity(tsVector, WHICH_CLASSIFY, island, (island+1));
	      vector<float> PhenDiv = getPhenotypeDiversity(tsVector, WHICH_CLASSIFY, island, (island+1));
	      fprintf(logFile2, "%f\t%f\t%i\t%i\t%f\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
		      i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3]);
	      //fprintf(logFile, "Most Fit: %f Average Fitness: %f of generation %i on island %i Standard Deviation: %f Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f\n",
	      //    i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	    }
	}
      if( (gen+1) % WHEN_PRINT_DATA == 0 ) {
    	  vector<float> HamDiv = getHammingDiversityForAllNodes();
    	  vector<float> FitDiv = getFitnessDiversity(tsVector, WHICH_CLASSIFY, 0, NUM_ISLANDS);
    	  vector<float> PhenDiv = getPhenotypeDiversity(tsVector, WHICH_CLASSIFY, 0, NUM_ISLANDS);  // we use the training set untel it's time to use the full test set
    	  fprintf(logFile1, "\n%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
    	  	    		  mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3]);
    	  //fprintf(logFile, "<---- Most Fit: %f on island %i at generation %i. Overall Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f. Phenotype: %f Fitness Div: %f ---->\n",
	    	//	  mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3], PhenotypeDiv, FitnessDiv[3]);
      }
      else if ((WHEN_PRINT_DATA > WHEN_FULL_TEST) && (((gen+1) % WHEN_FULL_TEST) == 0)) {
    	  fprintf(logFile1, "\n\t\t\t\t\t\t\t\t\t\t\t\t"); // make sure the full test sets line up
      }

      if( (gen+1) > WHEN_PRINT_DATA * 100) WHEN_PRINT_DATA *= 10;
      
      if (((gen+1) % WHEN_FULL_TEST) == 0) // need the +1 because population's gen-count has been updated during doOneGen
	{
	  mostFit = 0;
	  for(int island=0; island < NUM_ISLANDS; island++)
	    {
	      i = islands[island];
	      i.a_pop->updatePopulationFitness(all_tests,WHICH_CLASSIFY);
	      if(i.a_pop->getPopulationMaxFitness() > mostFit)
		{
		  mostFit = i.a_pop->getPopulationMaxFitness();
		  mostFitIsland = island;
		}
	      vector<float> HamDiv = i.a_pop->getInternalHammingDiversity();
	      vector<float> FitDiv = getFitnessDiversity(all_tests, WHICH_CLASSIFY, island, (island+1));
	      vector<float> PhenDiv = getPhenotypeDiversity(all_tests, WHICH_CLASSIFY, island, (island+1));
	      fprintf(logFile2, "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t%f\t%f\t%i\t%i\t%f\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
	      	    	  i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3]);
	      //fprintf(logFile2, "Full Testset: Most Fit: %f Average Fitness: %f of generation %i on island %i Standard Deviation: %f Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f\n",
	    	//  i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	      //i.a_pop->getBestIndividual().dumpConfMat(logFile2);
	    }
	   vector<float> HamDiv = getHammingDiversityForAllNodes();
 	   vector<float> FitDiv = getFitnessDiversity(all_tests, WHICH_CLASSIFY, 0, NUM_ISLANDS);
 	   vector<float> PhenDiv = getPhenotypeDiversity(all_tests, WHICH_CLASSIFY, 0, NUM_ISLANDS);
	   fprintf(logFile1, "\t\t%f\t%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
		   mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) HamDiv[0], (int) HamDiv[1], HamDiv[2], HamDiv[3], FitDiv[0], FitDiv[1], FitDiv[2], FitDiv[3], PhenDiv[0], PhenDiv[1], PhenDiv[2], PhenDiv[3]);
	   //fprintf(logFile1, "<---- Most Fit EVALUATED OVER ALL TESTS: %f on island %i at generation %i. Overall Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f. ---->",
	   //		   mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	}
      if ( gen+1 > WHEN_FULL_TEST * 10 ) WHEN_FULL_TEST *= 10;
    }

  // print end time
  end = time(NULL);
  fprintf(logFile1, "\n\nElapsed time %f seconds or %f minutes", difftime(end, start), difftime(end, start)/60.0 );
  fprintf(logFile2, "\n\nElapsed time %f seconds or %f minutes", difftime(end, start), difftime(end, start)/60.0 );
  fclose(logFile1);
  fclose(logFile2);
  return 0;
}

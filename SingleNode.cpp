#include "testinstance2.h"
#include "testset.h"
#include "individual2.h"
#include "population.h"
#include "SingleNode.h"
#include <time.h>
#include <stdlib.h>


using namespace std;

// GLOBAL VARIABLES
FILE *logFile;
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
      // can come back later and add the best, worst etc. if we want
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

/*
 * Computes the variance of the fitness of the individuals in the entire population.
 * Assumes that the fitness has been updated across the popualtions.
 */
float getFitnessVariance(vector<TestInstance2> testInstances, char classification) {
	int counter = 0;
	float mean, x, delta, var = 0.0;
	for (int i = 0; i < NUM_ISLANDS; ++i) {
		islands[i].a_pop->updatePopulationFitness(testInstances, classification);
		for (int j = 0; j < POP_SIZE; ++j) {
			++counter;
			x = islands[i].getIndividual(j).getFitness();
			delta = x - mean;
			mean = mean + (delta / (float)counter);
			var = var + (delta * (x - mean));
		}
	}
	var = var / (float)counter;
	return var;
	}

/**
 * Computes the fitness of each individual within each island across the given test set.
 * Returns a three-dimensional vector of [islands][individuals][test instances]
 */
vector< vector< vector<float> > > getPhenotypeDiversity(vector<TestInstance2> testInstances, char classification) {
	// use the full test set
	// then for the furst island
	// then for the first individual
	// compute how that individuals performs on each of the test cases
	// create a Vector that records how it did on each test case
	// repeat for all individuals and islands, adding a new nested vector
	// return that
		SingleNode tempIsland;
		Individual2 tempIndividual;
		vector< vector< vector<float> > > diversetableH (NUM_ISLANDS, vector< vector<float> >(POP_SIZE, vector<float>(testInstances.size()))); // for hi-fi classification
		//vector< vector< vector<int> > > diversetableL (NUM_ISLANDS, vector< vector<int> >(POP_SIZE, vector<int>(testInstances.size())));  // for low-fi classification, but we need to return something else
		for (int i = 0; i < NUM_ISLANDS; ++i) {
			tempIsland = islands[i];
			for (int j = 0; j < POP_SIZE; ++j) {
				tempIndividual = tempIsland.getIndividual(j);
				tempIndividual.resetConfMat();
				for (int k = 0; k < testInstances.size(); ++k) {
					switch (classification) {
						case 'h': diversetableH[i][j][k] = tempIndividual.classiHiFi(testInstances[k]); break;
						case 'l': diversetableH[i][j][k] = (float)tempIndividual.classify(testInstances[k]); break;
					}
				}
			}
		}
		return diversetableH;
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

  //initialize the logFile
  string s, outFile;
  stringstream out;

  bool depthdivision;
  if (argc == 15)
    {depthdivision = true;}
  else 
    {depthdivision = false;}
  
  outFile += "log/sim.";

  out << "rep-" << argv[2] << WHICH_FITNESS << "fifit-" << WHICH_CLASSIFY << "ficlsfy-" \
      << WHICH_MIGRATION << "ms" \
      << WHEN_MIGRATE << "mi" \
      << NUM_MIGRANTS_PER_ISLAND << "mn" \
      << NUM_ISLANDS << "i" \
      << NUM_NEIGHBORS << "n" \
      << POP_SIZE << "p" \
      << NUM_TEST_CASES_TO_USE << "ti" \
      << MAX_GENERATIONS << "g";
  s = out.str();
  outFile += s;

  remove(outFile.c_str());
  logFile = fopen(outFile.c_str(), "w");


  //run the generations
  //printf("Starting run with %d neighbors, %d migrants per island, %c migration strategy\n", NUM_NEIGHBORS, NUM_MIGRANTS_PER_ISLAND, WHICH_MIGRATION);
  SingleNode i;
  float mostFit=0;
  int mostFitIsland=-1;

  cout << "Counting generations: ";
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
	      vector<float> diversity = i.a_pop->getInternalHammingDiversity();
	      fprintf(logFile, "Most Fit: %f Average Fitness: %f of generation %i on island %i Standard Deviation: %f Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f\n",
		      i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	    }
	}
      if( (gen+1) % WHEN_PRINT_DATA == 0 ) {
    	  vector<float> diversity = getHammingDiversityForAllNodes();
	      fprintf(logFile, "<---- Most Fit: %f on island %i at generation %i. Overall Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f. ---->\n",
	    		  mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
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
	      vector<float> diversity = i.a_pop->getInternalHammingDiversity();
	      fprintf(logFile, "Full Testset: Most Fit: %f Average Fitness: %f of generation %i on island %i Standard Deviation: %f Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f\n",
	    	  i.a_pop->getPopulationMaxFitness(), i.a_pop->getPopulationAvgFitness(), i.a_pop->getGeneration(), island, i.a_pop->getStdev(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	      i.a_pop->getBestIndividual().dumpConfMat(logFile);
	    }
	   vector<float> diversity = getHammingDiversityForAllNodes();
	   //vector< vector< vector<float> > > phenotype = getPhenotypeDiversity(all_tests, WHICH_CLASSIFY);
	   //float performance = getFitnessVariance();
	   fprintf(logFile, "<---- Most Fit EVALUATED OVER ALL TESTS: %f on island %i at generation %i. Overall Max Diversity: %i Min Diversity: %i Average Diversity: %f Diversity Variance: %f. ---->\n",
		   mostFit, mostFitIsland, i.a_pop->getGeneration(), (int) diversity[0], (int) diversity[1], diversity[2], diversity[3]);
	}
      if ( gen+1 > WHEN_FULL_TEST * 10 ) WHEN_FULL_TEST *= 10;
    }

  // print end time
  end = time(NULL);
  fprintf(logFile, "Elapsed time %f seconds or %f minutes\n", difftime(end, start), difftime(end, start)/60.0 );
  fclose(logFile);
  return 0;
}

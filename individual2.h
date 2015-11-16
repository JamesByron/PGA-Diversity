#include <string>
#include <vector>
#include <time.h>
#include "config.h"

//typedef unsigned char byte;
using namespace std;

/**
   Individual2: revised to use binary comparison for efficiency

   Instances of this class represent individuals from a population.
   Stores an individual's genetic code in 'rule' and most recently computed
   fitness in 'fitness'.  
 */
class Individual2
{
public:
  Individual2();
  Individual2(string str);
  Individual2(int i);
  // classify a single test case
  //float fitnessHiFi(TestInstance2* ti);
  //float classiHiFi(TestInstance2* ti);
  //int classify(TestInstance2* ti);
  // updateFitnees with respect to all test cases
  //void updateFitness(TestSet* ts, char WHICH_FITNESS);
  //  void updateFitness(TestInstance2 * v);
  //void updateFitness(vector<TestInstance2>* v, char WHICH_CLASSIFY);
  //  void updateFitness(TestInstance2 * v, int rank, int numnodes);
  //void updateFitness(TestInstance2 testsplit[], int ntests, char WHICH_CLASSIFY);
  float getFitness() { return fitness; }
  // Testing
  //void findAccuracy(TestSet * ts, char WHICH_CLASSIFY);
  float getAccuracy() { return accuracy; }
  void breedNCross(Individual2 * kids, Individual2 ind);
  // single-point crossover reproduction
  void breed1Cross(Individual2 * kids, Individual2 ind);
  // two-point crossover reproduction
  void breed2Cross(Individual2 * kids, Individual2 ind);
  void setRule(string s);
  void setRule(unsigned char * ucs);
  void setRandomRule();
  unsigned char * getRule() { return rule; }
  unsigned int * getIntRule() { return intRule; }
  string getStringRule();
  void mutate();
  bool isSelected() { return selected;}
  void select() { selected = true; }
  void unselect() { selected = false;}
  void dumpConfMat(FILE *lf);
  void sortNums(int * cpts, int j);
  void resetConfMat();
  void updateDiversityRelevance(float relevance);
  float getDiversityRelevance() {return myDiversityRelevance; }
  void resetIntRule();

  void countFeats(signed char * featcounts, TestInstance2 * ti);
  float confMat[RULE_CASES][RULE_CASES];
  unsigned char rule[RULE_CASES*NUM_FEATURES];
  float fitness;

private:
  void auxBreedNCross(Individual2 * kids, Individual2 ind, int crossThisTime);
  unsigned char toUChar(string s);
  void splitbytes(unsigned char * n1, unsigned char * n2, unsigned char r1, unsigned char r2, int split);
  string byteToString(unsigned char c);
  //void resetConfMat(); moved to public
  //void sortNums(int * cpts, int j);
  bool selected;
  float accuracy;
  unsigned int intRule[RULE_CASES*NUM_FEATURES*8]; //  = {0};
  float myDiversityRelevance;
};

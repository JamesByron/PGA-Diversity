#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include "config.h"

using namespace std;

class Population
{
public:
  Population(int n, int nprocs, TestSet ts, int popsize);
  // accessors
  int getGeneration() { return generation; }
  float getPopulationMaxFitness() { return maxFitness; }
  float getPopulationAvgFitness() { return totalFitness/POP_SIZE; }
  float getStdev() {return stdev;}
  Individual2 getBestIndividual() {return bestIndiv;}
  // update fitness of all individuals
  void updatePopulationFitness(char WHICH_FITNESS);
  void updatePopulationFitness(vector<TestInstance2> allti, char WHICH_CLASSIFY);
  void populationAccuracy(char WHICH_CLASSIFY);
  void updatePopulationRelavance(vector<float>* relavance);
  // select individuals for . . . .
  void selectToSurvive(int n);
  void selectRandToMigrate(Individual2 * migrants, int num_migrants);
  void selectStrongToMigrate(Individual2 * migrants, int num_migrants);
  void selectWeakToMigrate(Individual2 * migrants, int num_migrants);
  // add immigrating individuals
  void processImmigrants(Individual2 * v, int n);
  void processOneImmigrant(Individual2 i);
  // generate offspring
  void generateOffspring(int n);
  // switch to next generation
  void nextGeneration(int PROB_MUTATE);
  Individual2* getIndividual(int index);
  vector<float> getInternalHammingDiversity();
  //vector<float> getExternalPopulationDiversity(Individual2 input);
  vector<int> calculateHammingForAll(Individual2* input);
  void addPopulationBitTotal(float * totals);
  void updatePopulationIntRules();
private:
  // functions
  int calculateHammingPair(unsigned int * a, unsigned int * b);
  //void updateInternalHammingDist();
  void unselectAll();
  int selectIndividual(int availablepop);
  int tournamentSelect(int availablepop);
  int altSelectIndividual();
  int selectWeakIndividual();
  int relavanceTournamentSelect(int availablepop);
  // variables
  int POP_SIZE;
  //vector<vector<int> > HAMMING_DIST;
  TestSet myTestSet;
  int generation;
  float avgFitness;
  float totalFitness;
  float totalInverseFitness;
  float maxFitness;
  float stdev;
  int myrank;
  int numprocs;
  int newpop_count;
  int genomeLength = NUM_FEATURES*RULE_CASES*8;
  Individual2 bestIndiv;
  //Individual2 b1[POP_SIZE];
  //Individual2 b2[POP_SIZE];
  Individual2 * b1;
  Individual2 * b2;
  Individual2 * mypop;
  Individual2 * newpop;
  //vector<float> relavanceVec;
};

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include "config.h"
#include "dataset.h"

#ifndef POPULATION_H
#define POPULATION_H

using namespace std;

class Population
{
public:
  Population(int n, int nprocs, DataSet<KRKTestInstance> ts, int popsize);
  // accessors
  int getGeneration() { return generation; }
  float getPopulationMaxFitness() { return maxFitness; }
  float getPopulationAvgFitness() { return totalFitness/POP_SIZE; }
  float getStdev() {return stdev;}
  Individual2 getBestIndividual() {return bestIndiv;}
  // update fitness of all individuals
  void updatePopulationFitness(char WHICH_FITNESS, char whichSet);
  void updateFitness(DataSet<KRKTestInstance>* ts, Individual2* individual, char WHICH_FITNESS, char whichSet);
  void updatePopulationRelevance(vector<float>* relevance);
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
  vector<int> calculateHammingForAll(Individual2* input);
  void addPopulationBitTotal(float * totals);
  void updatePopulationIntRules();
  DataSet<KRKTestInstance> myDataSet;
private:
  // functions
  int calculateHammingPair(unsigned int * a, unsigned int * b);
  void unselectAll();
  int selectIndividual(int availablepop);
  int tournamentSelect(int availablepop);
  int altSelectIndividual();
  int selectWeakIndividual();
  int relevanceTournamentSelect(int availablepop);
  // variables
  int POP_SIZE;
  int generation;
  float avgFitness;
  float totalFitness;
  float totalInverseFitness;
  float maxFitness;
  float stdev;
  int myrank;
  int numprocs;
  int newpop_count;
  // use RULE_LEN instead: int genomeLength = NUM_FEATURES*RULE_CASES*8;
  Individual2 bestIndiv;
  Individual2 * b1;
  Individual2 * b2;
  Individual2 * mypop;
  Individual2 * newpop;
};

#endif

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

class SingleNode
{
 public:
  SingleNode(int r, TestSet ts);
  SingleNode();
  //~Island() {};
  //void addToCustoms(vector<Individual2> v);
  //void addToCustoms(Individual2 * v, int n);
  void doOneGeneration(int thisgen);
  int sendMigrants();
  //  void NetFitnessAssessment(vector<TestInstance2> allti);
  // public variables
  Population * a_pop;
  Individual2 * customs;
  Individual2 getIndividual(int index);
  void addIslandBitTotal(float * totals);
  void updateNodeRelavance(vector<float>* islandRelavance);
  // int compareIndividualToPopulation(Individual2 input);
 private:
  bool selectIndividual(int ind);
  int myrank;
};

#include <string>
#include "config.h"
#include <vector>
#include "krktestinstance.h"

#ifndef DATASET_H
#define DATASET_H

using namespace std;

class DataSet
{
public:
  DataSet(vector<TestInstance*>* allti, int num, float split);
  DataSet() {}
  TestInstance * getTI(int i);
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(TestInstance ** ti, vector<TestInstance*> tests);
  void shuffle(vector<TestInstance*> * v);
  TestInstance ** test;
  TestInstance ** train;
};

#endif

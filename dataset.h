#include <string>
#include "config.h"
#include <vector>
#include "krktestinstance.h"

#ifndef DATASET_H
#define DATASET_H

using namespace std;

template <class T>
class DataSet
{
public:
  DataSet(vector<T*>* allti, int num, float split);
  DataSet() {}
  T * getTI(int i);
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(T ** ti, vector<T*> tests);
  void shuffle(vector<T*> * v);
  T ** test;
  T ** train;
};

#endif

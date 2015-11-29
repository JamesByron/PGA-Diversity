#include <string>
#include "config.h"
#include <vector>

using namespace std;

class DataSet
{
public:
  DataSet(vector<KRKTestInstance>* allti, int num, float split);
  DataSet() {}
  KRKTestInstance * getTI(int i)
  {
    if (i >= nTestInstances) {printf("getTestI: invalid index %d out of %d\n", i, nTestInstances); exit(-1);}
    return &test[i];
  }
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(KRKTestInstance * ti, vector<KRKTestInstance> tests);
  void shuffle(vector<KRKTestInstance> * v);
  KRKTestInstance * test;
};

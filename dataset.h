#include <string>
#include "config.h"
#include <vector>

using namespace std;

class DataSet
{
public:
  DataSet(vector<TestInstance2>* allti, int num, float split);
  DataSet() {}
  TestInstance2 * getTI(int i)
  {
    if (i >= nTestInstances) {printf("getTestI: invalid index %d out of %d\n", i, nTestInstances); exit(-1);}
    return &test[i];
  }
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(TestInstance2 * ti, vector<TestInstance2> tests);
  void shuffle(vector<TestInstance2> * v);
  TestInstance2 * test;
};

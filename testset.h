#include <string>
#include "config.h"
#include <vector>

using namespace std;

class TestSet
{
public:
  TestSet(vector<TestInstance2> allti, int n);
  TestSet(vector<TestInstance2> allti, int n, int depth);
  TestSet(int seed, vector<TestInstance2> allti, int n);
  TestSet(int seed, vector<TestInstance2> allti, int n, int depth);
  TestSet(vector<TestInstance2> allti, int num, float split);
  TestSet(int seed, vector<TestInstance2> allti, int num, float split);
  TestSet() {}
  //  ~TestSet()
  TestInstance2 getTI(int i)
  {
    if (i >= NUM_TEST_CASES_TO_USE) {printf("getTI: invalid index %d out of %d\n", i, NUM_TEST_CASES_TO_USE); exit(-1);}
    //return tset[i];
    return train[i];
  }
  TestInstance2 getTestI(int i)
  {
    if (i >= nTestInstances) {printf("getTestI: invalid index %d out of %d\n", i, nTestInstances); exit(-1);}
    return test[i];
  }
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  //float evaluateFitness(Individual2 ind);
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(TestInstance2 * ti, vector<TestInstance2> tests);
  void shuffle(vector<TestInstance2> v);
  //void  selectRandomTestInstances(vector<TestInstance2> tests);
  //vector<TestInstance2> fulltiset;
  TestInstance2 * tset;
  TestInstance2 * train;
  TestInstance2 * test;
};

#include "testinstance2.h"
//#include "individual2.h"
#include "testset.h"
#include <time.h>
#include <cmath>

using namespace std;

/** TestSet: a test-set instance should contain all of the testinstances that are known,
    but should provide distinct types of access to that total set:
    1. the full set
    2. all the instances of a particular classification
    3. a random sample of a given size from the full set
    4. a random sample of a given size from the instances of a particular classification
    5. a uniform sample across the different classes with the same number of each class

    Much of the work for providing such access should be done up front when creating the instance of TestSet.
    Need to preserve handle to the allti vector which gets passed in at creation, and then provide methods
    which get test instances from a particular one of the distinct samples.

    Perhaps, the TestSet class should provide these various testing services....

 */

TestSet::TestSet(vector<TestInstance2> allti, int num)
{
  //fulltiset = allti;
  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];
  // if (num == allti.size())
  // tset = allti;
  srand(time(NULL)+17033);
  shuffle(allti);
  //srand(93108);
  selectRandomTestInstances(tset, allti);
  //printf("sample has %d Testinsances\n", tset.size());
}

TestSet::TestSet(vector<TestInstance2> allti, int num, int depth)
{
  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];

  srand(time(NULL)+17033);
  //srand(93108);
  shuffle(allti);

  vector<TestInstance2> tempti;
  for(int i=0; i < allti.size(); i++)
    {
      if (allti[i].getDepth() == depth)
	tempti.push_back(allti[i]);
    }
  selectRandomTestInstances(tset, tempti);
  // ~vector (tempti);
}

TestSet::TestSet(int seed, vector<TestInstance2> allti, int num)
{
  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];
  // if (num == allti.size())
  // tset = allti;
  srand(seed);
  shuffle(allti);
  selectRandomTestInstances(tset, allti);
}

TestSet::TestSet(int seed, vector<TestInstance2> allti, int num, int depth)
{
  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];

  srand(seed);
  shuffle(allti);

  vector<TestInstance2> tempti;
  for(int i=0; i < allti.size(); i++)
    {
      if (allti[i].getDepth() == depth)
	tempti.push_back(allti[i]);
    }
  selectRandomTestInstances(tset, tempti);

  //~vector(tempti);
}

// constructors that creates two disjoint data sets for training and testing 
// with respective sizes given by the value of the split (which must correspond to num)

TestSet::TestSet(vector<TestInstance2> allti, int num, float split)
{
  int i = 0;
  if ( fabs(split - ((float) num)/allti.size()) > 0.1 ) {printf("BAD TRAIN/TEST SPLIT VALUE IN TestSet CONSTRUCTOR\n"); exit(-1);}

  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];
  // if (num == allti.size())
  // tset = allti;
  srand(time(NULL)+17033);
  shuffle(allti);
  //srand(93108);
  // get the training split
  train = new TestInstance2[num];
  for (int i=0; i < num; i++)
    train[i]=allti[i];
  // continue with the testing split
  test = new TestInstance2[allti.size() - num];
  for (int k=0; k < allti.size() - num; k++)
    test[k]=allti[i++];
  //printf("sample has %d Testinsances\n", tset.size());
}

TestSet::TestSet(int seed, vector<TestInstance2> allti, int num, float split)
{
  int i = 0;
  if ( fabs(split - ((float) num)/allti.size()) > 0.1 ) {printf("BAD TRAIN/TEST SPLIT VALUE IN TestSet CONSTRUCTOR\n"); exit(-1);}

  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti.size() - num;
  tset = new TestInstance2[num];
  // if (num == allti.size())
  // tset = allti;
  srand(seed);
  shuffle(allti);
  train = new TestInstance2[num];
  // get the training split
  for (i=0; i < num; i++)
    train[i]=allti[i];
  // continue with the testing split
  test = new TestInstance2[allti.size() - num];
  for (int k=0; k < allti.size() - num; k++)
    test[k]=allti[i++];
}


// SAMPLING FUNCTIONS

void TestSet::selectRandomTestInstances(TestInstance2 * ti, vector<TestInstance2> tests)
//void TestSet::selectRandomTestInstances(vector<TestInstance2> tests)
  // OUT ti: stuff the randomly selected cases into this space
  // IN tests: use these test instances to select a sample
{
  int index;
  
  for(int i=0; i<NUM_TEST_CASES_TO_USE; i++)
  {
    index = rand()%(tests.size());
    ti[i]=tests[index];
    //tset.push_back(tests[index]);
  }
}

void TestSet::shuffle(vector<TestInstance2> a)
{
  int r;
  TestInstance2 t;
  for (int n=0; n < 3; n++)
    for (int i=0; i < a.size(); i++)
      {
	r = rand() % a.size();
	t = a[r];
	a[r] = a[i];
	a[i] = t;
      }
}
      

// FITNESS EVALUATING

/*
float TestSet::evaluateFitness(Individual2 ind)
// Given an individual, determine its fitness
//  with respect to ... which test set?
//
{
  // for each testinstance ti
  //   call the individual's classification function on ti
  // set the individual's fitness via a new setter function in individual2
}
*/

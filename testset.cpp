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
// constructors that creates two disjoint data sets for training and testing
// with respective sizes given by the value of the split (which must correspond to num)
TestSet::TestSet(vector<TestInstance2> * allti, int num, float split)
{
  if ( fabs(split - ((float) num)/allti->size()) > 0.1 ) {printf("BAD TRAIN/TEST SPLIT VALUE IN TestSet CONSTRUCTOR\n"); exit(-1);}

  NUM_TEST_CASES_TO_USE = num;
  nTestInstances = allti->size();
  //nTestInstances = allti->size() - num;
  //train = new TestInstance2[num];
  test = new TestInstance2[allti->size()];
  srand(time(NULL)+17033);
  shuffle(allti);
  // get the training split
  int i = 0;
  while (i < allti->size()) {
	  //if (i < num) train[i] = (*allti)[i];
	  test[i] = (*allti)[i];
	  i++;
  }
  /*for (int i=0; i < num; i++)
    train[i]=(*allti)[i];
  // continue with the testing split
  test = new TestInstance2[allti->size() - num];
  for (int k=0; k < allti->size() - num; k++)
    test[k]=(*allti)[i++];*/
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
  }
}

void TestSet::shuffle(vector<TestInstance2> * a)
{
  int r;
  TestInstance2 t;
  for (int n=0; n < 3; n++)
    for (int i=0; i < a->size(); i++)
      {
	r = rand() % a->size();
	t = (*a)[r];
	(*a)[r] = (*a)[i];
	(*a)[i] = t;
      }
}

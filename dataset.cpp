//#include "testinstance.cpp"
#include "dataset.h"
#include <time.h>
#include <cmath>
#include <stdio.h>

using namespace std;

/** DataSet: a test-set instance should contain all of the krktestinstances that are known,
    but should provide distinct types of access to that total set:
    1. the full set
    2. all the instances of a particular classification
    3. a random sample of a given size from the full set
    4. a random sample of a given size from the instances of a particular classification
    5. a uniform sample across the different classes with the same number of each class

    Much of the work for providing such access should be done up front when creating the instance of DataSet.
    Need to preserve handle to the allti vector which gets passed in at creation, and then provide methods
    which get test instances from a particular one of the distinct samples.

    Perhaps, the DataSet class should provide these various testing services....

 */
// constructors that creates two disjoint data sets for training and testing
// with respective sizes given by the value of the split (which must correspond to num)
/*
template <class T>
DataSet<T>::DataSet(vector<T*> * allti, int num, float split) {
	//if ( fabs(split - ((float) num)/allti->size()) > 0.1 ) {printf("BAD TRAIN/TEST SPLIT VALUE IN DataSet CONSTRUCTOR\n"); exit(-1);}
	NUM_TEST_CASES_TO_USE = num;
	nTestInstances = allti->size() - num;
	train = new T*[num];
	test = new T*[allti->size() - num];
	srand(time(NULL)+17033);
	shuffle(allti);
	// get the training split
	int i = 0;
	while (i < num) {
		train[i] = (*allti)[i++];
	}
	// continue with the testing split
	for (int k=0; k < allti->size() - num; k++) {
		test[k]=(*allti)[i++];
	}
}
/*
template <class T>
T* DataSet<T>::getTI(int i) {
    if (i >= nTestInstances) {printf("getTestI: invalid index %d out of %d\n", i, nTestInstances); exit(-1);}
    return test[i];
  }
/*
// SAMPLING FUNCTIONS
template <class T>
void DataSet<T>::selectRandomTestInstances(T ** ti, vector<T*> tests)
//void DataSet::selectRandomTestInstances(vector<TestInstance> tests)
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
/*
template <class T>
void DataSet<T>::shuffle(vector<T*> * a)
{
	int r;
	T * t;
	for (int n=0; n < 3; n++)
		for (int i=0; i < a->size(); i++)
		{
			r = rand() % a->size();
			t = (*a)[r];
			(*a)[r] = (*a)[i];
			(*a)[i] = t;
		}
}*/

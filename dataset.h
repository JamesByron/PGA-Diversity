#include <string>
#include "config.h"
#include <vector>
#include <time.h>
#include <cmath>
#include <stdio.h>
#include "krktestinstance.h"

#ifndef DATASET_H
#define DATASET_H

using namespace std;

template <class T>
class DataSet
{
public:
  DataSet(vector<T*> * allti, int num, float split) {
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
  DataSet() {}
  T * getTI(int i) {
	    if (i >= nTestInstances) {printf("getTestI: invalid index %d out of %d\n", i, nTestInstances); exit(-1);}
	    return test[i];
	  }
  int trainSetSize() { return NUM_TEST_CASES_TO_USE; }
  int testSetSize() { return nTestInstances; }
  int NUM_TEST_CASES_TO_USE;
private:
  int nTestInstances;
  void selectRandomTestInstances(T ** ti, vector<T*> tests) {
		//void DataSet::selectRandomTestInstances(vector<TestInstance> tests)
		// OUT ti: stuff the randomly selected cases into this space
		// IN tests: use these test instances to select a sample
			int index;

			for(int i=0; i<NUM_TEST_CASES_TO_USE; i++)
			{
				index = rand()%(tests.size());
				ti[i]=tests[index];
			}
		}
  void shuffle(vector<T*> * a) {
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
		}
  T ** test;
  T ** train;
};

#endif

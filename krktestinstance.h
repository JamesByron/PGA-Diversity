#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "config.h"
//#include "dataset.h"
#include "testinstance.h"
//#include "individual2.h"
#ifndef KRKTESTINSTANCE_H
#define KRKTESTINSTANCE_H

using namespace std;

class KRKTestInstance: public TestInstance
{
public:
  KRKTestInstance() {};
  KRKTestInstance(string str);
  ~KRKTestInstance();
  string getStringRep();
  int classify(Individual2* individual);
  float fitnessHiFi(Individual2* individual);
  int distToNearCorner(int bkr, int bkf);
  int distToNearEdge(int rnk, int fl);
  void countFeats(signed char * featcounts, Individual2 * ind);
  string humanReadable() { return datasetForm; };
  unsigned char * getBinary() { return binary; };
  int getDepth() const { return depth; };
  int charDifference(char c, char d) { return (c-d)+1; };

protected:
	unsigned char binary[NUM_FEATURES];
	string datasetForm;
	int depth;
  int whiteKingFile;
  int whiteKingRank;
  int whiteRookFile;
  int whiteRookRank;
  int blackKingFile;
  int blackKingRank;
};

#endif

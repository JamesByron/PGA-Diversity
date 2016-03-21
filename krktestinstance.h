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

protected:
  int whiteKingFile;
  int whiteKingRank;
  int whiteRookFile;
  int whiteRookRank;
  int blackKingFile;
  int blackKingRank;
};

#endif

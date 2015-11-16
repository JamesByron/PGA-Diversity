#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "config.h"


using namespace std;

class TestInstance2
{
public:
  TestInstance2() {};
  TestInstance2(string str);
  ~TestInstance2();
  string humanReadable() { return datasetForm; };
  unsigned char * getBinary() { return binary; };
  int getDepth() const { return depth; };
  string getStringRep();

  // classify a single test case
  float fitnessHiFi(TestInstance2* ti);
  float classiHiFi(TestInstance2* ti);
  int classify(TestInstance2* ti);
  // updateFitnees with respect to all test cases
  void updateFitness(TestSet* ts, char WHICH_FITNESS);
  // void updateFitness(TestInstance2 * v);
  void updateFitness(vector<TestInstance2>* v, char WHICH_CLASSIFY);
  // void updateFitness(TestInstance2 * v, int rank, int numnodes);
  void updateFitness(TestInstance2 testsplit[], int ntests, char WHICH_CLASSIFY);
  void findAccuracy(TestSet * ts, char WHICH_CLASSIFY);

private:
  int charDifference(char c, char d) { return (c-d)+1; };
  string intToBinary(int i);
  string byteToString(unsigned char c);
  int distToNearCorner(int bkr, int bkf);
  int distToNearEdge(int rnk, int fl);
  int lowbyte(int i);
  int hibyte(int i);
  unsigned char binary[NUM_FEATURES];
  string datasetForm;
  int whiteKingFile;
  int whiteKingRank;
  int whiteRookFile;
  int whiteRookRank;
  int blackKingFile;
  int blackKingRank;
  int depth;
};

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "config.h"
#include "individual2.h"

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
  int classify(Individual2* individual);
  float fitnessHiFi(Individual2* individual);
private:
  int charDifference(char c, char d) { return (c-d)+1; };
  string intToBinary(int i);
  string byteToString(unsigned char c);
  int distToNearCorner(int bkr, int bkf);
  int distToNearEdge(int rnk, int fl);
  int lowbyte(int i);
  int hibyte(int i);
  void countFeats(signed char * featcounts, Individual2 * ind);
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

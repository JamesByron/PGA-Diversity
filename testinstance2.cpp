#include "testinstance2.h"
#include "dataset.h"
//#include "individual2.h"
#include <stdio.h>
#include <iostream>

using namespace std;

TestInstance2::TestInstance2(string str)
{
  datasetForm = str;
  whiteKingFile = charDifference(str[0], 'a');
  whiteKingRank = charDifference(str[2], '1');
  whiteRookFile = charDifference(str[4], 'a');
  whiteRookRank = charDifference(str[6], '1');
  blackKingFile = charDifference(str[8], 'a');
  blackKingRank = charDifference(str[10], '1');

  //set the depth

  switch(str[12]) {
  case 'd': depth = -1;
    //case 'd': depth = 16;
    break;
  case 'z': depth = 0;
    break;
  case 'o': depth = 1;
    break;
  case 't':
    if (str[13] == 'e') depth = 10;
    else if (str[14] == 'o') depth = 2;
    else if (str[14] == 'r') depth = 3;
    else if (str[14] == 'e') depth = 12;
    else depth = 13;
    break;
  case 'f':
    if (str[14] == 'v') depth = 5;
    else if (str[14] == 'f') depth = 15;
    else if (str.length() == 16) depth = 4;
    else depth = 14;
    break;
  case 's':
    if (str[13] == 'e') depth = 7;
    else if (str.length() == 15) depth = 6;
    else depth = 16;
    //else depth = -1;
    break;
  case 'e':
    if (str[13] == 'i') depth =8;
    else depth = 11;
    break;
  case 'n': depth = 9;
  }

  //set the binary
  // ranks and files 1..6
  string s = intToBinary(whiteKingFile);
  s += intToBinary(whiteKingRank);
  s += intToBinary(whiteRookFile);
  s += intToBinary(whiteRookRank);
  s += intToBinary(blackKingFile);
  s += intToBinary(blackKingRank);
  // absolute difference in rank 7..12
  /*
  s += intToBinary(1+abs(whiteKingRank - whiteRookRank));
  s += intToBinary(1+abs(whiteKingRank - blackKingRank));
  s += intToBinary(1+abs(whiteKingFile - whiteRookFile));
  s += intToBinary(1+abs(whiteKingFile - blackKingFile));
  s += intToBinary(1+abs(whiteRookRank - blackKingRank));
  s += intToBinary(1+abs(whiteRookFile - blackKingFile));
  */
  /*
  // augmented features 13..17
  int d2c = distToNearCorner(blackKingRank, blackKingFile);
  s += intToBinary(d2c);
  s += intToBinary(distToNearCorner(whiteKingRank, whiteKingFile));
  s += intToBinary(distToNearEdge(blackKingRank, blackKingFile));
  s += intToBinary(distToNearEdge(whiteKingRank, whiteKingFile));
  s += intToBinary((abs(whiteKingRank - blackKingRank)+abs(whiteKingFile - blackKingFile) > 8) ? 8 : abs(whiteKingRank - blackKingRank)+abs(whiteKingFile - blackKingFile));
  // features 18..21
  s += intToBinary(lowbyte(d2c+abs(whiteKingRank - blackKingRank)));
  s += intToBinary(hibyte(d2c+abs(whiteKingRank - blackKingRank)));
  s += intToBinary(lowbyte(d2c+abs(whiteKingFile - blackKingFile)));
  s += intToBinary(hibyte(d2c+abs(whiteKingFile - blackKingFile)));
  */
  /*
  printf("%d lo as: %s\n", (lowbyte(d2c+abs(whiteKingRank - blackKingRank))), intToBinary(lowbyte(d2c+abs(whiteKingRank - blackKingRank))).c_str());
  printf("%d hi as: %s\n", (hibyte(d2c+abs(whiteKingRank - blackKingRank))), intToBinary(hibyte(d2c+abs(whiteKingRank - blackKingRank))).c_str());
  printf("%d lo as: %s\n", (lowbyte(d2c+abs(whiteKingFile - blackKingFile))), intToBinary(lowbyte(d2c+abs(whiteKingFile - blackKingFile))).c_str());
  printf("%d hi as: %s\n", (hibyte(d2c+abs(whiteKingFile - blackKingFile))), intToBinary(hibyte(d2c+abs(whiteKingFile - blackKingFile))).c_str());
  */

  for(int f=0; f < NUM_FEATURES; f++)
    {
      unsigned char thischar = (unsigned char) 0;
      for(int i=0; i < 8; i++)
	{
	  thischar = thischar << 1;
	  if ( s[f*8+i] == '1' ) thischar++;
	}
      binary[f] = thischar;
    }

  return;
}

TestInstance2::~TestInstance2()
{
}

string TestInstance2::getStringRep()
{
  string s("");
  for (int i=0; i < NUM_FEATURES; i++)
    {
      s += byteToString(binary[i]);
      s += " ";
    }
  return s;
}

void TestInstance2::countFeats(signed char * featcounts, Individual2 * ind)
/** process the testinstance and stuff the feature match counts into the given array. Used by both fitnessHiFi and classiHiFi. */
{
	float result = 0.0;
	unsigned char * bin;
	for (int i=0; i < RULE_CASES; i++)
	{
		bin = getBinary();
		featcounts[i] = 0;
		for (int feat=0; feat < NUM_FEATURES; feat++)
			if( ((bin[feat] & ind->rule[i*NUM_FEATURES+feat]) != 0) || (!ind->rule[i*NUM_FEATURES+feat] && !bin[feat]) )
			{
				featcounts[i]++;
			}
	}
}

float TestInstance2::fitnessHiFi(Individual2* individual)
/** fine-grained fitness evaluation: looks at every feature of every classification rule (at a performance cost) accumulating
    the number of matched features for each rule.  Final fitness based on the distribution of values in some as-yet-undetermined manner.
*/
{
  int correctclass, correctmatched, bettermatched, samematched, mostmatched, nummostmatched;
  float result = 0.0;
  unsigned char * bin;
  signed char featcounts[RULE_CASES];

  countFeats(featcounts, individual);

  // now do something with the featcounts array
  correctclass = getDepth()+1;
  correctmatched = featcounts[correctclass];
  bettermatched = 0; samematched = 0; mostmatched = 0; nummostmatched = 0;
  for(int i=0; i < RULE_CASES; i++)
    {
      if ( featcounts[i] > mostmatched ) { mostmatched = featcounts[i]; nummostmatched = 1; }
      else if ( featcounts[i] == mostmatched ) nummostmatched++;
      if ( featcounts[i] > correctmatched ) bettermatched++;
      else if ( featcounts[i] == correctmatched ) samematched++;
    }
  // update the confusion matrix
  for (int i=0; i < RULE_CASES; i++)
    if ( featcounts[i] == mostmatched) individual->confMat[correctclass][i] += 1.0/nummostmatched;
  result = (((float) correctmatched)/NUM_FEATURES) * (1.0 / samematched) * (1.0 / (1 + bettermatched));
  //printf("fitnessHiFi: matched total of %d features with %d on correct rule, resulting in score: %f\n", matchedfeats, correctmatched, result);
  return result;
}


int TestInstance2::classify(Individual2* individual)
  /** Classify a single test instance as the first rule-case that completely matches all the features.
   */
{
  //printf("Entering classify\n");
  int correctclass = getDepth()+1;
  unsigned char * bin = getBinary();

  bool match;
  int i;
  for (i = 0; i < RULE_CASES; i++)
    {
      match = true;
      for (int j = 0; j < NUM_FEATURES; j++)
  {
    if ( ((bin[j] & individual->rule[i*NUM_FEATURES+j]) == 0) && (individual->rule[i*NUM_FEATURES+j] || bin[j]) )
      {
        match = false;
        break;
      }
  }
    if (match)
    { /* printf("Leaving classify\n");*/
      individual->confMat[correctclass][i]++; return (correctclass == i-1);} //classification
    }
  // Catch-all case
  i = rand() % RULE_CASES;
  individual->confMat[correctclass][i]++;
  //printf("Leaving classify\n");
  return (correctclass == i-1); 
}

string TestInstance2::intToBinary(int i)
// i is in the range [1,8]
// return an 8-digit string of 0s with a 1 in the i'th position from the right
{
  string s = "";
  while (s.length() < (8-i))
  { s = s + "0"; }
  if (s.length() < 8) 
    s = s + "1";
  while (s.length() < 8)
  { s = s + "0"; }
 
  return s;
}

string TestInstance2::byteToString(unsigned char c)
{
  string s("");
  int mask = 128;
  for (int i=0; i < 8; i++)
    {
      if ( (c & mask) == 0 )
	s += "0";
      else
	s += "1";
      mask /= 2;
    }
  return s;
}

int TestInstance2::distToNearCorner(int bkr, int bkf)
/** find the manhatten distance to the nearest corner of this rank and file
 */
{
  return min(bkr,9-bkr)+min(bkf,9-bkf);
}

int TestInstance2::distToNearEdge(int rnk, int fl)
/** find the distance to the nearest edge -- either along rank or file
 */
{
  return min(min(rnk,9-rnk),min(fl,9-fl));
}

int TestInstance2::lowbyte(int i)
/** for number between 1 and 16, return the number if less than or equal to 8, or 0 otherwise
 */
{
  if (i <= 8)
    return i;
  else
    return 0;
}

int TestInstance2::hibyte(int i)
/** for number between 1 and 16, return the amount above 8 if greater than 8, or 0 otherwise
 */
{
  if (i > 8)
    return i-8;
  else
    return 0;
}

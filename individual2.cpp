#include "individual2.h"
#include <stdio.h>
#include <iostream>

using namespace std;

Individual2::Individual2()
{
  selected = false;
}

Individual2::Individual2(int rl)
{
  string s("");
  if ( rl != RULE_LEN )  printf("w-ERROR: incorrect length of %d to Individual2(int) constructor\n", rl);
  for (int i = 0; i < RULE_LEN; i++)
    {
      if ( rand()%ONE_FREQ )
	s = s + "1";
      else
	s = s + "0";
    }
  setRule(s);
  selected = false;
}

Individual2::Individual2(string str)
  // expects a string of 0's and 1's, RULE_LEN chars long
{
  if (str.length() != RULE_LEN) printf("w-ERROR: incorrect length of %d to Individual2(string) constructor\n", (int)str.length());
  setRule(str);
  selected = false;
}

void Individual2::updateDiversityRelevance(float relevance) {
	myDiversityRelevance = relevance;
}

void Individual2::breedNCross(Individual2 * kids, Individual2 ind)
{
  int crossThisTime = (rand()%MAX_NUM_CROSS_BREED)+1;
  auxBreedNCross(kids, ind, crossThisTime);
}

void Individual2::breed1Cross(Individual2 * kids, Individual2 ind) { auxBreedNCross(kids, ind, 1); }

void Individual2::breed2Cross(Individual2 * kids, Individual2 ind) { auxBreedNCross(kids, ind, 2); }


void Individual2::setRule(string s)
{
  //printf("Entering setRule(string s)\n");
  for (int rc=0; rc < RULE_CASES; rc++)
    for (int f=0; f < NUM_FEATURES; f++)
	rule[rc*NUM_FEATURES+f] = toUChar(s.substr(rc*NUM_FEATURES*8+f*8, 8));
  //printf("Leaving setRule(sring s)\n");
}

void Individual2::setRule(unsigned char * ucs)
{
  //printf("Entering  setRule(unsigned char * ucs)\n");
  for (int i=0; i < RULE_CASES*NUM_FEATURES; i++)
    rule[i] = ucs[i];
  //printf("Leaving  setRule(unsigned char * ucs)\n");
}

void Individual2::setRandomRule()
{
  string s("");
  //printf("Entering setRule(string s)\n");
  for (int rc=0; rc < RULE_CASES; rc++)
    for (int f=0; f < NUM_FEATURES; f++)
      {
	s = "";
	for (int rdig=0; rdig < 8; rdig++)
	  if (rand()%ONE_FREQ)
	    s = s + "1";
	  else
	    s = s + "0";
	//rule[rc*NUM_FEATURES+f] = toUChar(s.substr(rc*NUM_FEATURES*8+f*8, 8));
	rule[rc*NUM_FEATURES+f] = toUChar(s);
      }
  //printf("Leaving setRule(sring s)\n");
}

string Individual2::getStringRule()
{
	exit(-1);
  string s("");
  for (int i=0; i < RULE_CASES*NUM_FEATURES; i++) {
    s += byteToString(rule[i]);
  }
  // check that it's right
  resetIntRule();
  for (int x = 0; x < s.size(); ++x) {
	  if ((s[x] == '1') && (1 != intRule[x])) {
		  cout << "there is a problem at " << x << " where s " << s[x] << " and rule " << intRule[x] << endl; cout.flush(); exit(-1);
	  }
	  else if ((s[x] == '0') && (0 != intRule[x])) {
		  cout << "there is a problem at " << x << " where s " << s[x] << " and rule " << intRule[x] << endl; cout.flush(); exit(-1);
	  }
  }
  cout << "good!" << endl;
  return s;
}

void Individual2::mutate()
  /** Revised 3/7/2009 (wfi): 
      Call this function ONLY IF the individual gets "hit", then roll the dice
      for which bit got hit, flip that bit, and that's it.
      NOTE: this precludes multiple bits getting flipped on the same generation.
   */
{
  //printf("Entering mutate\n");
  int hitPos = rand()%RULE_LEN;
  rule[hitPos/8] ^= ((unsigned char)1) << (unsigned char)(hitPos%8);
}

// PRIVATE

// auxBreedNCross ....
void Individual2::auxBreedNCross(Individual2 * kids, Individual2 ind, int crossThisTime)
{
  //printf("Entering auxBreedNCross\n");
  // initialize crossover points
  int cpts[MAX_NUM_CROSS_BREED];
  for (int k=0; k < crossThisTime; k++) cpts[k] = (rand()%(RULE_LEN-1))+1;
  sortNums(cpts, crossThisTime);
  //printf(". auxBNC: initialized array of %d crossover points:", crossThisTime);
  //for (int j=0; j < crossThisTime; j++) printf("%d ", cpts[j]); printf("\n");

  // continue with the new children
  unsigned char c1, c2, nc1 = 0, nc2 = 0, childRule1[RULE_CASES*NUM_FEATURES], childRule2[RULE_CASES*NUM_FEATURES];
  unsigned char * tmpcr, * cr1 = childRule1, * cr2 = childRule2;
  int n = 0, i=0;

  // for each cross-over point
  while(n < crossThisTime)
    {      
      // copy any bytes PRIOR to the byte containing the nth cross-over point
      while( ((i+1) * 8) < cpts[n] )
	{
	  cr1[i] = rule[i];
	  cr2[i] = ind.getRule()[i];
	  i++;
	}

      //split on this crossover point
      if ( (n+1 == crossThisTime && cpts[n]/8 == cpts[n-1]/8 ) // if this is last and prev was on this same byte
	   || ( n > 0 && cpts[n]/8 == cpts[n-1]/8 )) // or if this is not the first and the prev was on this same byte
	{
	  i--;
	  c1 = cr2[i];  // note: the rules cr1 and cr2 were swapped on the previous time
	  c2 = cr1[i];
	}
      else
	{
	  c1 = rule[i];
	  c2 = ind.getRule()[i];
	}
      splitbytes(&nc1, &nc2, c1, c2, cpts[n] % 8);
      cr1[i] = nc1; cr2[i] = nc2;
      nc1=0; nc2=0;

      // switch rules and move to the next cross point
      tmpcr = cr1;
      cr1 = cr2;
      cr2 = tmpcr;
      n++; i++;
    }
  //printf(". auxBNC: finished while loop -- copy remaining bytes\n");
  // copy any remaining bytes AFTER the byte containing the last of the crossThisTime crossover points
  while ( (i * 8) < RULE_LEN )
    {
      cr1[i] = rule[i];
      cr2[i] = ind.getRule()[i];
      i++;
    }

  //printf(". auxBNC: set the kids to the childRules\n");
  kids[0].setRule(childRule1);
  kids[1].setRule(childRule2);
  //printf("Leaving auxBreedNCross\n");
}


unsigned char Individual2::toUChar(string s)
{
  unsigned char auc = 0;
  for (int i=0; i < 8; i++)
    {
      auc = auc << 1;
      if ( s[i] == '1' ) auc++;
    }
  return auc;
}

void Individual2::splitbytes(unsigned char * n1, unsigned char * n2, unsigned char r1, unsigned char r2, int split)
  // IN-OUT n1 and n2 expected to be 0 at start
  // stuff n1 with the first "split" bits from r1, and n2 with the first "split" bits from r2
  // and then the rest of the bits of n1 from the rest of r2, and n2 from r1 respectively
{
  //printf("Entering splitbytes\n");
  unsigned char mask = 128;
  unsigned char tmp;
  //printf("BEFORE SPLIT on %d:: n1: %d, r1: %s, n2: %d, r2: %s\n",split, *n1, byteToString(r1).c_str(), *n2, byteToString(r2).c_str());
  for (int i=0; i < 8; i++)
    {
      if ( i == split ) { tmp = r2; r2 = r1; r1 = tmp; }
      *n1 = *n1 | (r1 & mask);
      *n2 = *n2 | (r2 & mask);
      mask >>= 1;
    }
  //printf("Leaving splitbytes\n");
  //printf("AFTER SPLIT      :: n1: %s, r1: %d, n2: %s, r2: %d\n",byteToString(*n1).c_str(), r1, byteToString(*n2).c_str(), r2);
}

string Individual2::byteToString(unsigned char c)
{
  string s("");
  unsigned char mask = 128;
  for (int i=0; i < 8; i++)
    {
      if ( !(c & mask) )
	s += "0";
      else
	s += "1";
      mask >>= 1;
    }
  return s;
}

/**
 * At each generation, resetIntRule is called after completing recombination and selection for the population.
 * When this is called, intRule is updated with the integer representations of the buts in the rule.  The intRule therefore
 * contains egight times as many elements as the rule itself, but the values in the int rule mirror the bits in the rule so that we
 * can add these values together to compute hamming distance.
 * Possible future work includes eliminating the intRule and moving the bit shifting to the population.cpp wherever intRule is used.
 */
void Individual2::resetIntRule() {
	static unsigned char mask;
	for (int i = 0; i < NUM_FEATURES*RULE_CASES; ++i) {
		mask = 128;
		for (int j = 0; j < 8; ++j) {
			if (rule[i] & mask) {
				intRule[(i*8)+j] = 1;
			}
			else intRule[(i*8)+j] = 0;
			mask >>= 1;
		}
	}
}

void Individual2::resetConfMat()
  // reset the confusion matrix
{
  for(int i=0; i < RULE_CASES; i++)
    for(int j=0; j < RULE_CASES; j++)
      confMat[i][j]=0;
}
    
void Individual2::dumpConfMat(FILE *lf)
  /** Print out this individual's confusion matrix as most recently populated
   */
{
  for(int i=0; i < RULE_CASES ; i++)
    {
      fprintf(lf, "ConfMat %3d: %8.2f", i, confMat[i][0]);
      for(int j=1; j < RULE_CASES; j++)
	fprintf(lf, ", %8.2f", confMat[i][j]);
      fprintf(lf, "\n");
    }
}


void Individual2::sortNums(int * cpts, int j)
// sort the j crossover points 
{
  while (j > 0){
    for (int i=0; i < (j-1); i++)
      {
	if (cpts[i] > cpts[j-1])
	  {
	      int temp = cpts[i];
	      cpts[i] = cpts[j-1];
	      cpts[j-1] = temp;
	  }
      }
    j--;
  }
}

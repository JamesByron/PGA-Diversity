//#include "dataset.h"
//#include <stdio.h>
//#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <string>
//#include "config.h"
#include "individual2.h"

using namespace std;

class TestInstance {
public:
	// these are the virtual classes that need to be implementd from this interface
	virtual string getStringRep() = 0;
	virtual void countFeats(signed char * featcounts, Individual2 * ind) = 0;
	virtual float fitnessHiFi(Individual2* individual) = 0;
	virtual int classify(Individual2* individual) = 0;

	// these are the pre-implemented functions that are standard across all testinstance implementations.
	string intToBinary(int i) {
		// i is in the range [1,8]
		// return an 8-digit string of 0s with a 1 in the i'th position from the right
		string s = "";
		while (s.length() < (8-i))
		{ s = s + "0"; }
		if (s.length() < 8)
			s = s + "1";
		while (s.length() < 8)
		{ s = s + "0"; }

		return s;
	}

	string byteToString(unsigned char c) {
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

	int lowbyte(int i) {
		/** for number between 1 and 16, return the number if less than or equal to 8, or 0 otherwise */
		if (i <= 8)
			return i;
		else
			return 0;
	}

	int hibyte(int i) {
		/** for number between 1 and 16, return the amount above 8 if greater than 8, or 0 otherwise */
		if (i > 8)
			return i-8;
		else
			return 0;
	}

	string humanReadable() { return datasetForm; };
	unsigned char * getBinary() { return binary; };
	int getDepth() const { return depth; };
	int charDifference(char c, char d) { return (c-d)+1; };

protected:
	unsigned char binary[NUM_FEATURES];
	string datasetForm;
	int depth;
};

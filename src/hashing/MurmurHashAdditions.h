/*
 * MurmurhasAdditions.h
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#ifndef MURMURHASHADDITIONS_H_
#define MURMURHASHADDITIONS_H_

#include "MurmurHash3.h"

#include "Snap.h"

const int MURMURSEED = 65765745;

const char base16[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };

inline TStr myhash(TStr in) {
	char* cstr = in.CStr();
	int len = in.Len();
	unsigned char* val = (unsigned char*) malloc(16);
	MurmurHash3_x86_128(cstr, len, MURMURSEED, val);
	char* stringVal = (char*) malloc(32 + 1);
	for (int i = 0; i < 16; ++i) {
		unsigned char c = val[i];
		unsigned char c1_index = c >> 4;
		unsigned char c2_index = c & 15;
		char c1 = base16[c1_index];
		char c2 = base16[c2_index];
		stringVal[2 * i] = c1;
		stringVal[2 * i + 1] = c2;
	}
	stringVal[32] = 0;
	return TStr(stringVal);
}

#endif /* MURMURHASHADDITIONS_H_ */

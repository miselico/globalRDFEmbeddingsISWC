/*
 * addPredicateTranspose.cpp
 *
 *  Created on: Mar 13, 2017
 *      Author: cochez
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

namespace {
typedef double real;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
} CREC;

int numberOfNodes = 7451907;

int firstPredicateIndex = numberOfNodes; //== the total number of nodes considered because they wil get IDs [0-amount]

void addPredicateTranspose() {

	FILE *fin = stdin;
	FILE *fout = stdout;

	CREC record;

	unsigned int i = 0;

	while (1) {
		int count = fread(&record, sizeof(CREC), 1, fin);
		if ((count == 0) && (feof(fin) != 0)) {
			break;
		} else if (count == 0) {
			std::cerr << "Error reading input, exiting" << std::endl;
			exit(2);
		}
		i++;
		if (record.word2 >= firstPredicateIndex) {
			//make extra record
			CREC extra { record.word2, record.word1, record.val };
			fwrite(&extra, sizeof(CREC), 1, fout);
		}

		fwrite(&record, sizeof(CREC), 1, fout);
	}
	fflush(fout);

}

} //end anonymous namespace

int NOTmain(int argc, char **argv) {
	addPredicateTranspose();
	return 0;
}

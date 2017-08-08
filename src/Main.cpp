/*
 * Main.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "Snap.h"

#include "RandomWalkExperiments.h"

#include "RDF2Co_occurence.h"

#include <iostream>

using namespace std;


int main(int argc, char **argv) {
	try{
	RDF2CO::performExperiments();
	} catch (char const* str){
		cout << str << endl;
		throw str;
	}
	return 0;
}

int mainRandomWalk(int argc, char **argv) {





	//RandomWalkExperiments::performExperiments(1);

	//return 0;






	if (argc < 2) {
		cerr << "Usage: commandName [number], where number from" << endl;
		cerr << "	1. Uniform" << endl << "	2. Predicate freq" << endl << "	3. Predicate inv. freq" << endl << "	4. Object freq" << endl << "	5. Object inv. freq" << endl << "	6. P-O freq" << endl
				<< "	7. P-O inv. freq" << endl << "	8. Pagerank weight" << endl << "	9. inv Pagerank weight" << endl << "	10. Pagerank split weight" << endl << "	11. inv Pagerank split weight" << endl
				<< "	12. Object inv. freq split" << endl << "	Object freq split == uniform"
				<< endl;
		return 1;
	}

	TStr strategy(argv[1]);
	int strategyNumber = strategy.GetInt();

	char * outFileName = NULL;
	if (argc > 2){
		outFileName = argv[2];
	}

	RandomWalkExperiments::performExperiments(strategyNumber, outFileName);



	//RDF2CO::performExperiments();
	return 0;
}


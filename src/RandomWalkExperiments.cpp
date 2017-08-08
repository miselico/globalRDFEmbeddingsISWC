/*
 * RandomWalkExperiments.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWalker.h"
#include "GraphWeigher.h"
#include "RandomWalkExperiments.h"
#include "RDF2Walk.h"
#include <iostream>

//FIXME remove include
#include "nTripleParser.h"

#include "Snap.h"

namespace RandomWalkExperiments {

using namespace std;

namespace {
static TStr RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static TStr OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");
}

void walkAroundStartFromStarts(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> >& graph, int walksPerNode, GraphWalker & walker, TextFileSink & sink, TVec<TInt> & starts) {
	long sumWalkLengths = 0;
	long generatedPaths = 0;
	int conceptCount = 0;
	long int nodeCount = graph->GetNodes();
	for (TVec<TInt>::TIter iter = starts.BegI(); iter < starts.EndI(); iter++) {
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI startNode = graph->GetNI(iter->Val);
		conceptCount++;
		//Generate paths
		for (int j = 0; j < walksPerNode; ++j) {
			TVec<TStr> path = walker.performWalk(graph, startNode);
			sumWalkLengths += path.Len();
			sink.consume(path);
			generatedPaths++;
		}

	}
	cout << "Found " << conceptCount << " concepts out of " << nodeCount << " entities. Generated " << walksPerNode << " per concept, i.e., " << generatedPaths << " paths" << endl;
	cout << "All paths generated, Average path length (total length including start node and predicates) = " << double(sumWalkLengths) / double(generatedPaths) << endl;
}

void walkAroundStartFromOWLThings(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> >& graph, int walksPerNode, GraphWalker & walker, TextFileSink & sink) {
	long sumWalkLengths = 0;
	long generatedPaths = 0;
	int conceptCount = 0;
	long int nodeCount = graph->GetNodes();
	for (int i = 0; i < nodeCount; i++) {
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI candidateNode = graph->GetNI(i);
		for (int outEdgeNr = 0; outEdgeNr < candidateNode.GetOutDeg(); ++outEdgeNr) {
			TStr predicate = candidateNode.GetOutEDat(outEdgeNr).P();
			if (predicate == RDF_TYPE) {
				TStr object = candidateNode.GetOutNDat(outEdgeNr);
				if (object == OWL_THING) {
					conceptCount++;
					//Generate paths
					for (int j = 0; j < walksPerNode; ++j) {
						TVec<TStr> path = walker.performWalk(graph, candidateNode);
						sumWalkLengths += path.Len();
						sink.consume(path);
						generatedPaths++;
					}
					//paths generated for this node
					break;
				}
			}
		}
	}
	cout << "Found " << conceptCount << " concepts out of " << nodeCount << " entities. Generated " << walksPerNode << " per concept, i.e., " << generatedPaths << " paths" << endl;
	cout << "All paths generated, Average path length (total length including start node and predicates) = " << double(sumWalkLengths) / double(generatedPaths) << endl;
}

THash<TStr, TFlt> readDBPediaPageRanks(TStr tsvFile) {
	THash<TStr, TFlt> ranks;

	PSIn FInPt = TFIn::New(tsvFile);
	TStr line;

	while (FInPt->GetNextLn(line)) {
		if (line.IsWs()) {
			continue;
		}
		if (line.SearchCh('#', 0) == 0) {
			//comment
			continue;
		}
		TStr resource = TStr("<") + line.LeftOf('\t') + TStr(">");
		TFlt rank = line.RightOf('\t').GetFlt();

		ranks.AddDat(resource, rank);
	}

	ranks.Pack();

	return ranks;
}

int performExperiments(int strategyNumber, char* outFileName) {

	//general walk settings
	int length = 4;
	int walksPerNode = 250;
	int seed = 45645;
	TStr walksOutFileName;
	if (outFileName == NULL){
		walksOutFileName = TStr("walks_strategy_") + TInt::GetStr(strategyNumber) + TStr(".txt");
	} else{
		walksOutFileName = TStr(outFileName);
	}
	//inputfile settings
	TStr ntriplesFileName = "allData.nt";
	//TStr ntriplesFileName = "SmallTest8_multiplePO.nt";
	TStr pageRankFileName = "pagerank.tsv";

	GraphWeigher * weigher = nullptr;

	switch (strategyNumber) {
	case 1:
		weigher = new UniformWeigher();
		break;
	case 2:
		weigher = new PredicateFrequencyWeigher();
		break;
	case 3:
		weigher = new InversePredicateFrequencyWeigher();
		break;
	case 4:
		weigher = new ObjectFrequencyWeigher();
		break;
	case 5:
		weigher = new InverseObjectFrequencyWeigher();
		break;
	case 6:
		weigher = new PredicateObjectFrequencyWeigher();
		break;
	case 7:
		weigher = new InversePredicateObjectFrequencyWeigher();
		break;
	case 8: {	//scoping pr
		THash<TStr, TFlt> pr = readDBPediaPageRanks(pageRankFileName);
		weigher = new PushDownWeigher(pr, 0.2);
	}
		break;
	case 9: {	//scoping pr
		THash<TStr, TFlt> pr = readDBPediaPageRanks(pageRankFileName);
		for (THash<TStr, TFlt>::TIter iter = pr.BegI(); iter < pr.EndI(); iter++) {
			iter.GetDat() = 1.0 / iter.GetDat();
		}
		weigher = new PushDownWeigher(pr, 0.2);
	}
		break;
	case 10: {	//scoping pr
		THash<TStr, TFlt> pr = readDBPediaPageRanks(pageRankFileName);
		weigher = new SplitDownWeigher(pr, 0.2);
	}
		break;
	case 11: {	//scoping pr
		THash<TStr, TFlt> pr = readDBPediaPageRanks(pageRankFileName);
		for (THash<TStr, TFlt>::TIter iter = pr.BegI(); iter < pr.EndI(); iter++) {
			iter.GetDat() = 1.0 / iter.GetDat();
		}
		weigher = new SplitDownWeigher(pr, 0.2);
	}
		break;
	case 12:
		weigher = new InverseObjectFrequencyWeigherSplitDown();
		break;
	default:
		cerr << "Only numbers between 1 and 12 are valid";
		return 1;
	}

	if (weigher == nullptr) {
		cerr << "Strategy not implemented yet";
		return 1;
	}

	FILE* walksOutFile = fopen(walksOutFileName.CStr(), "w");
	TextFileSink sink(walksOutFile);
	cout << "Reading in data and weighing" << endl;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraphIgnoringLiterals(ntriplesFileName, *weigher);
//alternative: also leave out leafs in the graph -> also change a couple of lines below!!!

//	TPair<TPt<TNodeEdgeNet<TStr, WeightedPredicate> >, TVec<TInt> > grapAndStarts = buildWalkableGraphIgnoringLiteralsAndLeafs(ntriplesFileName, *weigher);

	delete weigher;

	cout << "Done weighing. Starting walks." << endl;

	RandomProportionalWalker walker(length, seed);


	walkAroundStartFromOWLThings(graph, walksPerNode, walker, sink);

//alternaive: also leave out leafs in the graph -> also change a couple of lines above!!!
//	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = grapAndStarts.Val1;
//	TVec<TInt> startingpoints = grapAndStarts.Val2;
//	walkAroundStartFromStarts(graph, walksPerNode, walker, sink, startingpoints);

	if (fclose(walksOutFile) == 0) {
		cout << "file closed correctly" << endl;
	} else {
		cerr << "Error closing file" << endl;
		exit(2);
	}

	return 0;
}

namespace {
//OLD experiments, not in use any longer

void experiment0() {
	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
	//TStr file = "sample-wikidata-terms-fragment.nt";
	//TStr file = "sample-wikidata-terms.nt";
	//TStr file = "SmallTest3.nt";
	TStr file = "SmallTest5_duplicates.nt";

	//char* outfile = "frequencies_output.bin";

	//	computeFrequencies(file, inverseFrequencyWeigher, fopen(outfile, "w"));

	const char* walksoutfile = "walks";

	FILE* f = fopen(walksoutfile, "w");

	int length = 8;
	int seed = 45645;
	RandomProportionalWalker walker(length, seed);
	LengthEnforcingWalker walker2(walker, 2 * length + 1);
	InversePredicateFrequencyWeigher weigher;
	TextFileSink sink = TextFileSink(f);

	int amount = 1E6;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraph(file, weigher);

	TRnd anRnd(24332);

	for (long int i = 0; i < amount; ++i) {
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI start = graph->GetRndNI(anRnd);
		TVec<TStr> path = walker.performWalk(graph, start);
		sink.consume(path);
	}

	fclose(f);
}

void experiment1() {
//	* Let's start with 4 hops i.e., initial node + 4 predicates + 4 nodes.
//	* Shorter walks could be allowed, but I don't think that is going to happen in the DBpedia graph.
//	* All entities should be included (~5M). We can ignore the literals. Last time I used the following files from the DBpedia download server: instance_types, instance_types_transitive, categories, mappingbased_objects, Wikipedia_links and external_links. I think we can use the same subset again.
//	* let's try the uniform probability first. If that works, we can then use the inverse frequency of the predicate, and then the PageRank of the object.

	TStr file = "allData.nt";
//	TStr file = "SmallTest6_literals.nt";

	const char* walksoutfile = "walks.txt";
	FILE* f = fopen(walksoutfile, "w");
	int length = 4;
	int walksPerNode = 200;
	int seed = 45645;
	RandomProportionalWalker walker(length, seed);
	//LengthEnforcingWalker walker2 (walker, 2 * length + 1);
	UniformWeigher weigher;
	TextFileSink sink(f);

	cout << "Reading in data" << endl;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraphIgnoringLiterals(file, weigher);

	walkAroundStartFromOWLThings(graph, walksPerNode, walker, sink);
	if (fclose(f) == 0) {
		cout << "file closed correctly" << endl;
	} else {
		cerr << "Error closing file" << endl;
		exit(2);
	}
}

void experiment2() {
	cout << "Using experimental weights now!" << endl;
	TStr DBpediaFile = "pagerank_smallTest_for4.tsv";
	THash<TStr, TFlt> pr = readDBPediaPageRanks(DBpediaFile);

//	TStr file = "allData.nt";
	TStr file = "SmallTest7_owlThings.nt";

	const char* walksoutfile = "walksPRBiased.txt";
	FILE* f = fopen(walksoutfile, "w");
	int length = 4;
	int walksPerNode = 1000;
	int seed = 45645;
	RandomProportionalWalker walker(length, seed);
	//LengthEnforcingWalker walker2 (walker, 2 * length + 1);
	PushDownWeigher weigher(pr, 0.2);
	TextFileSink sink(f);

	cout << "Reading in data" << endl;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraphIgnoringLiterals(file, weigher);

	walkAroundStartFromOWLThings(graph, walksPerNode, walker, sink);
	if (fclose(f) == 0) {
		cout << "file closed correctly" << endl;
	} else {
		cerr << "Error closing file" << endl;
		exit(2);
	}
}

}

}

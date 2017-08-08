/*
 * RDF2Co_occurence.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include <iostream>
#include <fstream>
#include <string>
#include "Snap.h"
#include "WeightedPredicate.h"
#include "nTripleParser.h"
#include "BCA.h"
#include "GraphWeigher.h"
#include "utils.h"
#include "boost/dynamic_bitset.hpp"
#include "boost/graph/graph_traits.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "boost/lexical_cast.hpp"
//#include <boost/graph/topological_sort.hpp>
#include <utility>                   // for std::pair
#include "PrintTime.h"

namespace BCAOrder {

class my_bitset: public boost::dynamic_bitset<> {

private:
	boost::dynamic_bitset<>::size_type lastFound;
	boost::dynamic_bitset<>::size_type lastSetToTrue;

public:

	my_bitset() :
			lastFound(0), lastSetToTrue(0) {

	}

	boost::dynamic_bitset<>::size_type find_any() {
		if ((*this)[lastSetToTrue]) {
			return lastSetToTrue;
		}
		this->lastFound = this->find_next(lastFound);
		if (lastFound != this->npos) {
			return lastFound;
		} else {
			lastFound = find_first();
			return lastFound;
		}

	}

	dynamic_bitset& setTrueAndRecord(size_type n) {
		this->set(n, true);
		this->lastSetToTrue = n;
		return *this;
	}

};

//does the vector contain all numbers from 0 till all.Len()-1 ?
bool allNumbersIn(TVec<TInt> all) {
	TVec<TInt> copy(all);
	copy.Sort();
	for (int i = 0; i < all.Len(); i++) {
		if (copy[i] != i) {
			return false;
		}
	}
	return true;
}

class Node {
public:
	//TODO this also works with std::set instead of unordered. Check what is faster for large graphs.
	int ID;
	boost::unordered_set<Node*> inedgeSources;
	boost::unordered_set<Node*> outedgeDestination;

	Node(int ID) :
			ID(ID), inedgeSources(5), outedgeDestination(5) {			//Intentionally empty
	}

	int getInDeg() const {
		return inedgeSources.size();
	}
	int getOutDeg() const {
		return outedgeDestination.size();
	}

};

// directed graph (single directed edge between an ordered pair of nodes)
//This implementation is made specifically for determineBCVcomputeOrder
class MyGraph: public std::vector<Node> {
public:
	MyGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
		cout << currentTime() << "copying graph" << endl;
		const int totalNodes = baseGraph->GetNodes();
		this->reserve(totalNodes);
		for (int i = 0; i < totalNodes; i++) {
			this->emplace_back(i);
		}
		cout << currentTime() << "copying graph - nodes copied" << endl;

		//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
		//We also remove all self edges
		for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
			int src = EI.GetSrcNId();
			int dst = EI.GetDstNId();
			if (src != dst) {
				//boost::unordered_set<Node*>::value_type  == Node*
				//eclipse does not find these defnitions, gcc does...
				boost::unordered_set<Node*>::value_type sourceNode = &(*this)[src];
				boost::unordered_set<Node*>::value_type destinationNode = &(*this)[dst];
				std::pair<boost::unordered::iterator_detail::c_iterator<boost::unordered::detail::ptr_node<Node*> >, bool> a;

				a = sourceNode->outedgeDestination.insert(destinationNode);

				if (a.second) {
					destinationNode->inedgeSources.insert(sourceNode);
				}
			}
		}
		cout << currentTime() << "graph copied" << endl;
	}
};

///**
// * Attempts to determine the order in which most BCA computations can be reused.
// * This is using the heuristic that
// *
// * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
// * 2. if no such node exists, then the node with highest in degree is selected
// *
// * the input graph MUST have nodes indexed 0 till Nodes()-1
// *
// */
TVec<TInt> determineBCAcomputeOrder(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	MyGraph prunable_graph(baseGraph);
	const int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original

	TVec<TInt> finalOrder;
	finalOrder.Reserve(startSize);

	my_bitset zeroOutDegrees;
	zeroOutDegrees.resize(startSize, false);

	boost::dynamic_bitset<> todo;
	todo.resize(startSize, true);

	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin(); node != prunable_graph.cend(); node++) {
		if (node->getOutDeg() == 0) {
			zeroOutDegrees[node->ID] = true;
		}
	}

	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		Node& toBeRemoved = prunable_graph[n];

		for (boost::unordered_set<Node*>::const_iterator dependant = toBeRemoved.inedgeSources.cbegin(); dependant != toBeRemoved.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			boost::unordered_set<Node*>::value_type toBeRemovedAdress = &toBeRemoved;
			(*dependant)->outedgeDestination.erase(toBeRemovedAdress);
		}
		//finalize, delete node
		//the only things remaining are some bookkeeping:
#ifndef NDEBUG
		//probably not necessary anyhow
		toBeRemoved.inedgeSources.clear();
#endif
		todo[n] = false;
		finalOrder.Add(n);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}
	}

	cout << currentTime() << "After first fast phase, " << finalOrder.Len() << "/" << startSize << " nodes are done, starting iterative phase" << endl;

	//now the more general case including loops is handled
	TMaxPriorityQueue<TInt> highestInDegree;
	//set-up the  highestInDegree PQ,

	for (std::vector<Node>::const_iterator node = prunable_graph.cbegin(); node < prunable_graph.cend(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		if (todo[node->ID]) {
			highestInDegree.Insert(node->ID, node->getInDeg() + 1);
		}
	}

	//algo start

	while (finalOrder.Len() < startSize) {
		//break at least one loop
		TInt k = -1;
		do {
			IAssert(highestInDegree.Size() > 0);
			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
			k = highestInDegree.PopMax();
		} while (!todo[k]);		//this check is needed because the priorities for nodes removed in the loop below are not removed from the queueu, only marked as done
		//k will be removed add all nodes which will get a zero out degree to set
		Node& kNode = prunable_graph[k];

		for (boost::unordered_set<Node*>::const_iterator dependant = kNode.inedgeSources.cbegin(); dependant != kNode.inedgeSources.cend(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			(*dependant)->outedgeDestination.erase(&kNode);
		}

		//update the priorities of all nodes k is pointing to
		for (boost::unordered_set<Node*>::const_iterator dest = kNode.outedgeDestination.cbegin(); dest != kNode.outedgeDestination.cend(); dest++) {
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that dest.id != k
			highestInDegree.SetPriority((*dest)->ID, (*dest)->getInDeg() - 1 + 1);
			//start removing already
			(*dest)->inedgeSources.erase(&kNode);
		}

		//finalize, the removal is already done in the adjecancy lists
#ifndef NDEBUG
		//probably not necessary anyhow
		kNode.inedgeSources.clear();
		kNode.outedgeDestination.clear();
#endif
		todo[k] = false;
		finalOrder.Add(k);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}

		if (finalOrder.Len() == startSize) {
			break;
		}

		//remove as much as possible without having to break a loop
		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			Node& toBeRemoved = prunable_graph[n];

			for (boost::unordered_set<Node*>::const_iterator dependant = toBeRemoved.inedgeSources.cbegin(); dependant != toBeRemoved.inedgeSources.cend(); dependant++) {
				//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
				if ((*dependant)->getOutDeg() == 1) {
					zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
				}
				//start removing already
				(*dependant)->outedgeDestination.erase(&toBeRemoved);
			}
			//finalize, delete node
			//the only things remaining are some bookkeeping:
#ifndef NDEBUG
			//probably not necessary anyhow
			toBeRemoved.inedgeSources.clear();
#endif
			todo[n] = false;
			finalOrder.Add(n);
			if (finalOrder.Len() % infoFrequency == 0) {
				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
			}
		}

	}

	cout << currentTime() << "All done" << finalOrder.Len() << "/" << startSize << " done" << endl;

	IAssert(startSize == finalOrder.Len());
	IAssert(todo.find_first() == todo.npos);
	IAssert(allNumbersIn(finalOrder));
	return finalOrder;
}

void ComputeBCAOrder(TStr inputFilename, TStr outputFileName) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputFilename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	cout << currentTime() << " Now computing BCV compute order" << endl;
	TVec<TInt, int> order;
	//order = determineBCVcomputeOrderOptimized(graph);
	order = determineBCAcomputeOrder(graph);
	cout << currentTime() << "end determining BCV compute order, writing to file" << endl;

	ofstream myfile(outputFileName.CStr());
	for (TVec<TInt>::TIter iter = order.BegI(); iter < order.EndI(); iter++) {
		myfile << int(*iter) << '\n';
	}
	myfile.flush();
	if (!myfile.good()) {
		throw "WTF";
	}
	myfile.close();
	cout << currentTime() << "done writing" << endl;
}

TVec<TInt> readBCAOrder(TStr precomputedBCAOrderFile, int expectedCount) {
	TVec<TInt> result;
	result.Reserve(expectedCount);
	ifstream myfile(precomputedBCAOrderFile.CStr());
	int next;
	while (myfile >> next) {
		result.Add(next);
	}
	myfile.close();

	if (result.Len() != expectedCount) {
		throw "Size of the precomputed BCA order file is not as expected";
	}

	return result;
}

}

namespace co_occurence_computer {
/**
 * typedefs compatible with the input expected by glove
 */
typedef double real;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
} CREC;

/**
 * Get a table which assigns a unique index to each (predicate, object) pair.
 *
 * The index will be strictly greater as zero since that is what glove uses to skip an entry.
 */
static THash<TPair<TStr, TStr>, TInt> createPairedWordIndexTable(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TPair<TStr, TStr>, TInt> table;
	int counter = 1;
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) {

		int node_i_outdeg = NI.GetOutDeg();
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			TStr predicate = NI.GetOutEDat(outEdge);
			TStr object = NI.GetOutNDat(outEdge);
			TPair<TStr, TStr> pair = TPair<TStr, TStr>(predicate, object);
			if (!table.IsKey(pair)) {
				table.AddDat(pair, counter);
				counter++;
			}
		}
	}
	return table;
}

/*
 * The wordIndexTable is indexed from 0 as it should be, but we need IDs which start from 1. This function abstracts this away
 */

static int graphIDToGloveID(int graphID) {
	return graphID + 1;
}

/**
 * For each unique predicate in the graph add a unique ID it to the returned hash.
 *
 * The returned vector contains the strings in the same order as their IDs
 *
 * The used IDs range from graph.Nodes() till graph.Nodes+(numberOfUniquePredicates=returnValue.Len()) exclusive
 */
static TPair<TVec<TStr>, THash<TStr, int>> computePredicateIDs(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TStr, int> preds;
	TVec<TStr> labels;
	unsigned int currentID = graph->GetNodes();
	for (int i = 0; i < graph->GetEdges(); ++i) {
		TStr label = graph->GetEDat(i);
		if (!preds.IsKey(label)) {
			preds.AddDat(label, currentID);
			labels.Add(label);
			currentID++;
		}
	}
	return TPair<TVec<TStr>, THash<TStr, int>>(labels, preds);
}

bool isEntity(TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI const & candidateNode) {
	for (int outEdgeNr = 0; outEdgeNr < candidateNode.GetOutDeg(); ++outEdgeNr) {
		const WeightedPredicate data = candidateNode.GetOutEDat(outEdgeNr);
		TStr predicate = data.P();
		if (predicate == RDF_TYPE) {
			TStr object = candidateNode.GetOutNDat(outEdgeNr);
			if (object == OWL_THING) {
				return true;
			}
		}
	}
	return false;
}

class Co_occurenceComputer {
private:
	//currently preceded by _ to prevent name clashes in not yet adapted methods
	TVec<TInt> _order;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > _weightedGraph;
	//TPair<TVec<TStr>, THash<TStr, int>> predLabelsAndIDs;
	THash<TStr, int> _predGraphIDs;
	TVec<TStr> _predicateLabels;

public:

	Co_occurenceComputer(TStr inputGraphFileName, GraphWeigher& weighingStrategy) {
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputGraphFileName);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
		_order = BCAOrder::determineBCAcomputeOrder(graph);
		_weightedGraph = weighingStrategy.weigh(graph);
		TPair<TVec<TStr, int>, THash<TStr, int> > predLabelsAndIDs = computePredicateIDs(graph);
		_predGraphIDs = predLabelsAndIDs.Val2;
		_predicateLabels = predLabelsAndIDs.Val1;
	}

	Co_occurenceComputer(TStr inputGraphFileName, TStr precomputedBCAOrderFile, GraphWeigher& weighingStrategy) {
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputGraphFileName);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
		_order = BCAOrder::readBCAOrder(precomputedBCAOrderFile, graph->GetNodes());
		_weightedGraph = weighingStrategy.weigh(graph);
		TPair<TVec<TStr, int>, THash<TStr, int> > predLabelsAndIDs = computePredicateIDs(graph);
		_predGraphIDs = predLabelsAndIDs.Val2;
		_predicateLabels = predLabelsAndIDs.Val1;
	}

	/**
	 *
	 * Compute the BCA score for each pair in the graph under the given weighing strategy.
	 * Additionally, adds a score for each edge as well.
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 *
	 *
	 */
	void computeFrequenciesIncludingEdges(double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize, bool onlyEntities) {

		const int infoFrequency = _order.Len() > 1000000 ? 100000 : 10000;
		THash<TInt, BCV> bcvCache;

		int counter = 0;
		for (TVec<TInt>::TIter iter = _order.BegI(); iter < _order.EndI(); iter++) {
			const int focusWordGraphID = iter->Val;
			//		//only take specific one:
			//		if (candidateNode.GetDat() != "<http://dbpedia.org/ontology/Province>"){
			//			continue;
			//		}

			if (onlyEntities) {
				if (!isEntity(_weightedGraph->GetNI(focusWordGraphID))) {
					continue;
				}
			}
			BCV combinedbcv = computeBCAIncludingEdgesCached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, _predGraphIDs, bcvCache);
			const int focusWordGloveID = graphIDToGloveID(focusWordGraphID);
			if (normalize) {
				combinedbcv.removeEntry(focusWordGraphID);
				combinedbcv.normalizeInPlace();
			}
			for (THash<TInt, TFlt>::TIter iter = combinedbcv.BegI(); iter < combinedbcv.EndI(); iter++) {
				int contextWordGraphID = iter.GetKey();
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
				//CREC crec = CREC { word1:contextWordGloveID, word2:focusWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}
			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/" << _order.Len() << " BCV computations" << endl;
			}
		}

		//still need to write all node labels to the vocab file
		for (int i = 0; i < this->_weightedGraph->GetNodes(); i++) {
			TStr nodeLabel = _weightedGraph->GetNDat(i);
			fprintf(glove_vocab_file_out, "%s nofr\n", nodeLabel.CStr());
		}
		//still need to write all predicates to the vocab file

		for (TVec<TStr>::TIter it = _predicateLabels.BegI(); it < _predicateLabels.EndI(); it++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", it->CStr());
		}
	}

	/**
	 *
	 * Compute the BCA score for each noe pair in the graph under the given weighing strategy. Ignores predicates.
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 */
	void computeFrequencies(double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize, bool onlyEntities) {
		cerr << "computeFrequencies is not yet tested thoroughly, check results with care";
		const int infoFrequency = _order.Len() > 1000000 ? 100000 : 10000;
		THash<TInt, BCV> bcvCache;
		int counter = 0;

		for (TVec<TInt>::TIter iter = _order.BegI(); iter < _order.EndI(); iter++) {
			const int focusWordGraphID = iter->Val;

			if (onlyEntities) {
				if (!isEntity(_weightedGraph->GetNI(focusWordGraphID))) {
					continue;
				}
			}
			BCV bcv = computeBCACached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, bcvCache);

			int focusWordGloveID = graphIDToGloveID(focusWordGraphID);
			if (normalize) {
				bcv.removeEntry(focusWordGraphID);
				bcv.normalizeInPlace();
			}

			for (THash<TInt, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
				int contextWordGraphID = iter.GetKey();
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}
			counter++;
			if ((counter % infoFrequency) == 0) {
				cout << currentTime() << "Processed " << counter << "/" << _order.Len() << " BCV computations" << endl;
			}
		}
		//still need to write all node labels to the vocab file
		for (int i = 0; i < this->_weightedGraph->GetNodes(); i++) {
			TStr nodeLabel = _weightedGraph->GetNDat(i);
			fprintf(glove_vocab_file_out, "%s nofr\n", nodeLabel.CStr());
		}
	}

};

class Co_occurenceComputer_Ultimate {
private:
	//currently preceded by _ to prevent name clashes in not yet adapted methods
	TVec<TInt> _order;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > _weightedGraph;
	//TPair<TVec<TStr>, THash<TStr, int>> predLabelsAndIDs;
	THash<TStr, int> _predGraphIDs;
	TVec<TStr> _predicateLabels;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > _weightedReverseGraph;
	TVec<TInt> _orderReverse;

public:

	Co_occurenceComputer_Ultimate(TStr inputGraphFileName, GraphWeigher& weighingStrategy, GraphWeigher & reverseWeighingStrategy) {
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputGraphFileName);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
		_order = BCAOrder::determineBCAcomputeOrder(graph);
		_weightedGraph = weighingStrategy.weigh(graph);
		TPair<TVec<TStr, int>, THash<TStr, int> > predLabelsAndIDs = computePredicateIDs(graph);
		_predGraphIDs = predLabelsAndIDs.Val2;
		_predicateLabels = predLabelsAndIDs.Val1;
		//reverse graph
		TPt<TNodeEdgeNet<TStr, TStr> > reversed = reverseGraph(graph);
		_orderReverse = BCAOrder::determineBCAcomputeOrder(reversed);
		_weightedReverseGraph = reverseWeighingStrategy.weigh(reversed);
	}

	/**
	 *
	 * Compute the BCA score for each pair in the graph under the given weighing strategy.
	 * Additionally, adds a score for each edge as well.
	 *
	 * Furthermore, it also performs a reverse walk and adds the result of that to the BCVs
	 * The reverse walk can be performed according to a different weighing strategy
	 *
	 * Outputs the score as a sparse matrix which can be fed to glove.
	 *
	 */
	void computeFrequenciesIncludingEdgesTheUltimate(double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize, bool onlyEntities) {
		const int infoFrequency = _order.Len() > 1000000 ? 100000 : 10000;
		THash<TInt, BCV> bcvForwardCache;
		{			//scoping forward
			int counter = 0;
			for (TVec<TInt>::TIter iter = _order.BegI(); iter < _order.EndI(); iter++) {
				const int focusWordGraphID = iter->Val;

				if (onlyEntities) {
					if (!isEntity(_weightedGraph->GetNI(focusWordGraphID))) {
						continue;
					}
				}
				//we only want the side effect of the BCV being added to the cache!
				computeBCAIncludingEdgesCached(_weightedGraph, focusWordGraphID, bca_alpha, bca_eps, _predGraphIDs, bcvForwardCache);
				counter++;
				if ((counter % infoFrequency) == 0) {
					cout << currentTime() << "Processed " << counter << "/" << _order.Len() << " BCV FORWARD computations " << endl;
				}
			}
		}
		{			//scoping backward
			int backwardCounter = 0;
			THash<TInt, BCV> bcvBackwardCache;
			for (TVec<TInt>::TIter iter = _orderReverse.BegI(); iter < _orderReverse.EndI(); iter++) {
				const int focusWordGraphID = iter->Val;

				if (onlyEntities) {
					if (!isEntity(_weightedReverseGraph->GetNI(focusWordGraphID))) {
						continue;
					}
				}
				BCV backwardBCV = computeBCACached(_weightedReverseGraph, focusWordGraphID, bca_alpha, bca_eps, bcvBackwardCache);

				//combine with what is in the forward cache
				BCV forwardBCV = bcvForwardCache.GetDat(focusWordGraphID);
				//remove from forward cache. This is not 100% necessary, but helps ensuring program correctness, and saves a bit of memory. Any forward entry must only be needed once.
				bcvForwardCache.DelKey(focusWordGraphID);

				forwardBCV.add(backwardBCV);
				BCV &combinedBCV = forwardBCV;

				if (normalize) {
					combinedBCV.removeEntry(focusWordGraphID);
					combinedBCV.normalizeInPlace();
				}
				const int focusWordGloveID = graphIDToGloveID(focusWordGraphID);

				for (THash<TInt, TFlt>::TIter iter = combinedBCV.BegI(); iter < combinedBCV.EndI(); iter++) {
					int contextWordGraphID = iter.GetKey();
					int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
					double freq = iter.GetDat();
					CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
					fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
				}
				backwardCounter++;
				if ((backwardCounter % infoFrequency) == 0) {
					cout << currentTime() << "Processed " << backwardCounter << "/" << _order.Len() << " BCV BACKWARD computations" << endl;
				}
			}
		}

		//still need to write all node labels to the vocab file
		for (int i = 0; i < this->_weightedGraph->GetNodes(); i++) {
			TStr nodeLabel = _weightedGraph->GetNDat(i);
			fprintf(glove_vocab_file_out, "%s nofr\n", nodeLabel.CStr());
		}
		//still need to write all predicates to the vocab file

		for (TVec<TStr>::TIter it = _predicateLabels.BegI(); it < _predicateLabels.EndI(); it++) {
			fprintf(glove_vocab_file_out, "%s nofr\n", it->CStr());
		}

	}
};

/////////////////////////not yet adapted (reusing parts + optiimzations) static methods////////////////////


/**
 * Implementation incomplete!!!
 */
void computeFrequenciesPushBCA(TStr filename, GraphWeigher& weighingStrategy, FILE *fout) {
	cerr << "computeFrequenciesPushBCA is not yet completely implemented";
	throw "not implemented";

	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraph(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

	THash<TStr, int> wordIndexTable = graphAndNodeIndex.Val2;

	THash<TPair<TStr, TStr>, TInt> pairwordIndexTable = createPairedWordIndexTable(graph);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);

	graph.Clr();

	cout << "done weighing" << endl;

	cerr << "TODO : check - does the indexing for glove have to start from 1 or 0 ??";
	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {
		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
			continue;
		}
		PBCV bcv = computePBCA(weightedGraph, i, 0.10, 0.000001);
		int subjectIndex = wordIndexTable.GetDat(weightedGraph->GetNDat(i));
		for (THash<TPair<TInt, TInt>, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
			WeightedPredicate wpred = weightedGraph->GetEDat(iter.GetKey().Val1);
			TStr pred = wpred.P();
			TStr obj = weightedGraph->GetNDat(iter.GetKey().Val2);
			int offset = wordIndexTable.Len();
			int pairIndex = pairwordIndexTable.GetDat(TPair<TStr, TStr>(pred, obj)) + offset;
			double freq = iter.GetDat();

			CREC crec = CREC { word1:subjectIndex, word2:pairIndex, val: freq };
			fwrite(&crec, sizeof(CREC), 1, fout);

			cout << subjectIndex << " has " << pairIndex << " freq " << freq << endl;
		}
	}

	//	TTmStopWatch w (true);
	//	int needed = 10000;
	//	int * selected = (int*) malloc(needed * sizeof(int));
	//	int skipped = 0;
	//	for (int i = 0; i < needed + skipped && i < weightedGraph->GetNodes(); ++i) {
	//		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
	//			++skipped;
	//			continue;
	//		} else {
	//			selected[i - skipped] = i;
	//		}
	//	}
	//	for (int index = 0; index < needed && index < weightedGraph->GetNodes(); index++) {
	//		int id = selected[index];
	//		PBCV bcv = computePBCA(weightedGraph, id, 0.10, 0.000000000001);
	//		//cout << bcv.toString(weightedGraph) << endl;
	//		if (index % 1000 == 0) {
	//			cout << "another 1000" << weightedGraph->GetNDat(id).CStr() << "->" << bcv.toString(weightedGraph) << endl;
	//		}
	//	}
	//
	//	w.Stop();
	//	cout << w.GetMSecInt() << "ms" << endl;
	return;
}

}		//end namespace co_occurence_computer

namespace RDF2CO {

//precompute the BCA order
//void performExperiments() {
//	TStr graphInputFile = "../../datasets/dbPedia/allData27_30M.nt";
//	auto BCAOrderFile = "BCAComputeOrder.txt";
//	BCAOrder::ComputeBCAOrder(graphInputFile, BCAOrderFile);
//}

void performExperiments() {
	TStr graphInputFile = "../../datasets/dbPedia/allData27_30M.nt";

//TStr graphInputFile = "SmallTest.nt";
//	TStr graphInputFile = "SmallTest9_loop.nt";

	UniformWeigher weigher;
//two options, if the BCA order has been precomputed:
//TStr precomputedBCAOrderFile = "BCAComputeOrder.txt";
//Co_occurenceComputer c(graphInputFile, precomputedBCAOrderFile, weigher);
//if it has not been precomputed:
	co_occurence_computer::Co_occurenceComputer c(graphInputFile, weigher);
//	co_occurence_computer::Co_occurenceComputer_Ultimate c(graphInputFile, weigher, weigher);
//now, c can be used to compute co_occurence matrices
	for (double alpha = 0.5; alpha <= 0.5; alpha += 0.1) {
		for (double eps = 0.001; eps >= 0.00001; eps /= 10) {

			string glove_input_file = "glove_input_file_out_alpha_" + boost::lexical_cast<std::string>(alpha) + "_eps_" + boost::lexical_cast<std::string>(eps) + ".bin";

			string glove_vocab_file = "glove_vocab_file_out_alpha_" + boost::lexical_cast<std::string>(alpha) + "_eps_" + boost::lexical_cast<std::string>(eps) + ".bin";

			cout << "writing to " << glove_input_file << endl;
			cout << "\tand " << glove_vocab_file << endl;

			FILE* glove_input_file_out = fopen(glove_input_file.c_str(), "w");
			FILE* glove_vocab_file_out = fopen(glove_vocab_file.c_str(), "w");

			//c.computeFrequenciesIncludingEdgesTheUltimate(alpha, eps, glove_input_file_out, glove_vocab_file_out, true, false);

			c.computeFrequenciesIncludingEdges(alpha, eps, glove_input_file_out, glove_vocab_file_out, true, false);

			fclose(glove_input_file_out);
			fclose(glove_vocab_file_out);
		}
	}

}

void performExperimentsOLD() {
//	TStr file = "wikidata-simple-statements-10_000000-sample.nt";
//TStr file = "sample-wikidata-terms-fragment.nt";
//TStr file = "sample-wikidata-terms.nt";
	TStr file = "../../datasets/dbPedia/allData27_30M.nt";
//TStr file = "SmallTest8_multiplePO.nt";

	FILE* glove_input_file_out = fopen("glove_input_file_out.bin", "w");
	FILE* glove_vocab_file_out = fopen("glove_vocab_file_out", "w");

	auto BCAOrderFile = "BCAComputeOrder.txt";
	BCAOrder::ComputeBCAOrder(file, BCAOrderFile);

//InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

//computeFrequenciesPushBCA(file, weigher, outfile);
//	bool normalize = true;
//	bool onlyEntities = false;
//	computeFrequenciesIncludingEdges(file, weigher, 0.80, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize, onlyEntities);
//computeFrequencies(file, weigher, 0.50, 0.05, glove_input_file_out, glove_vocab_file_out);

//computeFrequenciesIncludingEdgesTheUltimate(file, weigher, weigher, 0.70, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize);

	fclose(glove_input_file_out);
	fclose(glove_vocab_file_out);
}

}


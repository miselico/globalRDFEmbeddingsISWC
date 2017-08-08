/*
 * RDF2Co_occurence_slow.cpp
 *
 *  Created on: Mar 9, 2017
 *      Author: cochez
 */

#include <iostream>
#include "Snap.h"
#include "boost/dynamic_bitset.hpp"
#include "PrintTime.h"

//slow versions of algorithms, kept for later reference
namespace {

using namespace std;

/**
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 * First version without much of any optimization
 *
 */
TVec<TInt> determineBCVcomputeOrderOLD(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	TPt<TNEGraph> prunable = TNEGraph::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		prunable->AddEdge(EI.GetSrcNId(), EI.GetDstNId());
	}
	baseGraph = NULL; //making sure we do not accidentally write to the original
	TVec<TInt> result;

	while (prunable->GetNodes() > 0) {

		TVec<TInt> withZeroOutDegree;
		for (TNEGraph::TNodeI iter = prunable->BegNI(); iter < prunable->EndNI(); iter++) {
			if (iter.GetOutDeg() == 0) {
				withZeroOutDegree.Add(iter.GetId());
			}
		}
		if (withZeroOutDegree.Len() > 0) {
			//found some
			result.AddV(withZeroOutDegree);
			//remove from graph
			for (TVec<TInt>::TIter ID = withZeroOutDegree.BegI(); ID < withZeroOutDegree.EndI(); ID++) {
				prunable->DelNode(ID->Val);
			}
			continue;
		}
		//none with zero out degree. Take the one with highest indegree.
		int highestIndegree = -1;
		int withHighestInDegree = -1;
		for (TNEGraph::TNodeI iter = prunable->BegNI(); iter < prunable->EndNI(); iter++) {
			if (iter.GetInDeg() > highestIndegree) {
				highestIndegree = iter.GetInDeg();
				withHighestInDegree = iter.GetId();
			}
		}
		if (withHighestInDegree == -1) {
			cerr << "error, none found with in degree";
			exit(5);
		}
		result.Add(withHighestInDegree);
		prunable->DelNode(withHighestInDegree);
	}
	return result;
}


/**
 * Second version with some successful optimizations. The main bottleneck is copying the graph for large graphs.
 *
 *
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 */
TVec<TInt> determineBCVcomputeOrder(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
	//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
	//We also remove all self edges
	TPt<TNGraph> prunable = TNGraph::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		if (src != dst) {
			prunable->AddEdge(src, dst);
		}
	}
	int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original
	TVec<TInt> finalOrder;
	THashSet<TInt> zeroOutDegNodes;
	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		if (node.GetOutDeg() == 0) {
			zeroOutDegNodes.AddKey(node.GetId());
		}
	}
	while (zeroOutDegNodes.Len() != 0) {
		const TInt n = *(zeroOutDegNodes.BegI());
		zeroOutDegNodes.DelKey(n);
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		TNGraph::TNodeI niter = prunable->GetNI(n);
		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
			int mid = niter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI m = prunable->GetNI(mid);
			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
			if (m.GetOutDeg() == 1) {
				zeroOutDegNodes.AddKey(m.GetId());
			}
		}
		//finalize
		prunable->DelNode(n);
		finalOrder.Add(n);
	}
	//now the more general case including loops is handled

	TMaxPriorityQueue<TInt> highestInDegree;
	//set-up the  highestInDegree PQ,
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
	}
	//algo start
	while (prunable->GetNodes() > 0) {
		while (zeroOutDegNodes.Len() != 0) {
			const TInt n = *(zeroOutDegNodes.BegI());
			zeroOutDegNodes.DelKey(n);
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			TNGraph::TNodeI niter = prunable->GetNI(n);
			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
				int mid = niter.GetInNId(inEdgeNumber);
				TNGraph::TNodeI m = prunable->GetNI(mid);
				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
				if (m.GetOutDeg() == 1) {
					zeroOutDegNodes.AddKey(m.GetId());
				}
			}
			//set indegree of n to 0 in PQueueu. PQueue does not support direct removal
			highestInDegree.SetPriority(n, 0.0);
			//finalize
			prunable->DelNode(n);
			finalOrder.Add(n);
		}
		if (prunable->GetNodes() == 0) {
			break;
		}
		IAssert(highestInDegree.Size() > 0);
		//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
		IAssert(highestInDegree.GetMaxPriority() > 0.0);
		cerr.flush();
		TInt k = highestInDegree.PopMax();
		//add all nodes which will get a zero out degree to set
		IAssert(prunable->IsNode(k));
		TNGraph::TNodeI kiter = prunable->GetNI(k);
		IAssert(kiter.GetOutDeg() > 0);
		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
			int lid = kiter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			if (l.GetOutDeg() == 1) {
				zeroOutDegNodes.AddKey(l.GetId());
			}
		}
		//update the priorities of all nodes k is pointing to
		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
			int lid = kiter.GetOutNId(outEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
		}
		//finalize
		prunable->DelNode(k);
		finalOrder.Add(k);
	}
	IAssert(startSize == finalOrder.Len());
	return finalOrder;
}

class my_bitset: public boost::dynamic_bitset<> {

private:
	boost::dynamic_bitset<>::size_type lastFound = 0;
	boost::dynamic_bitset<>::size_type lastSetToTrue = 0;

public:
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


/**
 * Further optmimizations. using custom bitset.
 *
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 * the input graph MUST have nodes indexed 0 till Nodes()-1
 *
 */
TVec<TInt> determineBCVcomputeOrderOptimized(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
	//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
	//We also remove all self edges
	cout << currentTime() << "copying graph" << endl;
	TPt<TNGraph> prunable = TNGraph::New();
	prunable->Reserve(baseGraph->GetNodes(), baseGraph->GetEdges());
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	cout << currentTime() << "copying graph - nodes copied" << endl;
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		if (src != dst) {
			prunable->AddEdge(src, dst);
		}
	}
	const int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original
	cout << currentTime() << "graph copied" << endl;

	TVec<TInt> finalOrder;
	finalOrder.Reserve(startSize);

	my_bitset zeroOutDegrees;
	zeroOutDegrees.resize(startSize, false);

	boost::dynamic_bitset<> todo;
	todo.resize(startSize, true);

	//vector<bool> zeroOutDegrees;
	//zeroOutDegrees.

	//THashSet<TInt> zeroOutDegNodes;
	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		if (node.GetOutDeg() == 0) {
			zeroOutDegrees[node.GetId()] = true;
		}
	}
	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		//const int n = zeroOutDegrees.find_first();
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		TNGraph::TNodeI niter = prunable->GetNI(n);
		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
			int mid = niter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI m = prunable->GetNI(mid);
			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
			if (m.GetOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord(m.GetId());
				//zeroOutDegrees[m.GetId()] = true;
			}
		}
		//finalize
		prunable->DelNode(n);
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
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
	}
	//algo start
	while (prunable->GetNodes() > 0) {
		//while (zeroOutDegrees.any()) {
		//	int n = zeroOutDegrees.find_first();
		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			TNGraph::TNodeI niter = prunable->GetNI(n);
			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
				int mid = niter.GetInNId(inEdgeNumber);
				TNGraph::TNodeI m = prunable->GetNI(mid);
				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
				if (m.GetOutDeg() == 1) {
					zeroOutDegrees[m.GetId()] = true;
				}
			}
			//finalize
			prunable->DelNode(n);
			//We do not set indegree of n to 0 in PQueueu. PQueue does not support direct removal, but this is somewhat expensive. We use the to_do bitset instead.
			//			highestInDegree.SetPriority(n, 0.0);

			todo[n] = false;
			finalOrder.Add(n);
			if (finalOrder.Len() % infoFrequency == 0) {
				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
			}
		}
		if (prunable->GetNodes() == 0) {
			break;
		}
		TInt k = -1;
		do {
			IAssert(highestInDegree.Size() > 0);
			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
			k = highestInDegree.PopMax();
		} while (!todo[k]);
		//add all nodes which will get a zero out degree to set
		IAssert(prunable->IsNode(k));
		TNGraph::TNodeI kiter = prunable->GetNI(k);
		IAssert(kiter.GetOutDeg() > 0);
		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
			int lid = kiter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			if (l.GetOutDeg() == 1) {
				zeroOutDegrees[l.GetId()] = true;
			}
		}
		//update the priorities of all nodes k is pointing to
		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
			int lid = kiter.GetOutNId(outEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
		}
		//finalize
		prunable->DelNode(k);
		todo[k] = false;
		finalOrder.Add(k);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}

	}
	IAssert(startSize == finalOrder.Len());
	IAssert(todo.find_first() == todo.npos);
	IAssert(allNumbersIn(finalOrder));
	return finalOrder;
}



// * Uses boost graphs as the Snap ones seem to slow when they get large, never materialized as Michael Cochez could not get the boost graph library to work properly
//
///**
// * Some potential future ideas for optimization:
// * do we need to know about both in AND out edges or only about one of them if we only need to know the count of the in/out edges, we can store that in an array!
// */
//
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
//TVec<TInt> determineBCVcomputeOrderOptimizedBOOST(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
//	bla();
//
//	using namespace boost;
//// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
////We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
////We also remove all self edges
//
//	typedef adjacency_list<hash_setS, // Store out-edges of each vertex in a unordered set
//			listS,
//			// Store vertex set in a std::vector
//			bidirectionalS // The file dependency graph is directed
//	> Graph;
//
//	typedef graph_traits<Graph> GT;
//	typedef GT::vertex_descriptor Vertex;
//	typedef GT::edge_descriptor Edge;
//	typedef GT::vertex_iterator vertex_iter;
//
////Assumption: basegraph has nodes 0->nodes-1
//
//	cout << currentTime() << "copying graph" << endl;
//
//	Graph prunable_graph(baseGraph->GetNodes());
//
//	const unsigned int startSize = baseGraph->GetNodes();
//	{ //scoping allVertices
//	  //We keep iterators, they will even remain valid on removal, as we us a linked list!
//		std::vector<vertex_iter> allVertices(startSize);
//
//		int i = 0;
//		for (std::pair<vertex_iter, vertex_iter> vp = vertices(prunable_graph); vp.first != vp.second; ++vp.first) {
//			allVertices[i] = vp.first;
//			i++;
//		}
//		IAssert(allVertices.size() == startSize);
//		cout << currentTime() << "copying graph - nodes copied" << endl;
//		for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
//			int src = EI.GetSrcNId();
//			int dst = EI.GetDstNId();
//
//			Vertex srcV = *(allVertices[src]);
//			Vertex dstV = *(allVertices[dst]);
//
//			if (src != dst) {
//				add_edge(srcV, dstV, prunable_graph);
//				//prunable->AddEdge(src, dst);
//			}
//		}
//
//	} //end scoping allVertices
//	baseGraph = NULL; //making sure we do not accidentally write to the original
//	cout << currentTime() << "graph copied" << endl;
//
//	TVec<TInt> finalOrder;
//	finalOrder.Reserve(startSize);
//
//	my_bitset zeroOutDegrees;
//	zeroOutDegrees.resize(startSize, false);
//
//	boost::dynamic_bitset<> todo;
//	todo.resize(startSize, true);
//
//	//we do a first quicker check to eliminate all nodes which are not part of a cycle
//	typename property_map<Graph, vertex_index_t>::type index = get(vertex_index, prunable_graph);
//	for (std::pair<vertex_iter, vertex_iter> vp = vertices(prunable_graph); vp.first != vp.second; ++vp.first) {
//		unsigned int outDegree = out_degree(*vp.first, prunable_graph);
//
//		if (outDegree == 0) {
//			int a = index[*vp.first];
//
//			//zeroOutDegrees[]
//		}
//
//	}
//
////	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
////		if (node.GetOutDeg() == 0) {
////			zeroOutDegrees[node.GetId()] = true;
////		}
////	}
//
////	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;
////
////	unsigned long int n;
////	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
////		//const int n = zeroOutDegrees.find_first();
////		zeroOutDegrees[n] = false;
////		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
////		TNGraph::TNodeI niter = prunable->GetNI(n);
////		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
////			int mid = niter.GetInNId(inEdgeNumber);
////			TNGraph::TNodeI m = prunable->GetNI(mid);
////			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
////			if (m.GetOutDeg() == 1) {
////				zeroOutDegrees.setTrueAndRecord(m.GetId());
////				//zeroOutDegrees[m.GetId()] = true;
////			}
////		}
////		//finalize
////		prunable->DelNode(n);
////		todo[n] = false;
////		finalOrder.Add(n);
////		if (finalOrder.Len() % infoFrequency == 0) {
////			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////		}
////	}
////
////	cout << currentTime() << "After first fast phase, " << finalOrder.Len() << "/" << startSize << " nodes are done, starting iterative phase" << endl;
////
////	//now the more general case including loops is handled
////	TMaxPriorityQueue<TInt> highestInDegree;
////	//set-up the  highestInDegree PQ,
////	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
////		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
////		//We add one to indicate that even a node with 0 in degree is still a valid node.
////		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
////	}
////	//algo start
////	while (prunable->GetNodes() > 0) {
////		//while (zeroOutDegrees.any()) {
////		//	int n = zeroOutDegrees.find_first();
////		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
////			zeroOutDegrees[n] = false;
////			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
////			TNGraph::TNodeI niter = prunable->GetNI(n);
////			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
////				int mid = niter.GetInNId(inEdgeNumber);
////				TNGraph::TNodeI m = prunable->GetNI(mid);
////				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
////				if (m.GetOutDeg() == 1) {
////					zeroOutDegrees[m.GetId()] = true;
////				}
////			}
////			//finalize
////			prunable->DelNode(n);
////			//We do not set indegree of n to 0 in PQueueu. PQueue does not support direct removal, but this is somewhat expensive. We use the to_do bitset instead.
////			//			highestInDegree.SetPriority(n, 0.0);
////
////			todo[n] = false;
////			finalOrder.Add(n);
////			if (finalOrder.Len() % infoFrequency == 0) {
////				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////			}
////		}
////		if (prunable->GetNodes() == 0) {
////			break;
////		}
////		TInt k = -1;
////		do {
////			IAssert(highestInDegree.Size() > 0);
////			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
////			k = highestInDegree.PopMax();
////		} while (!todo[k]);
////		//add all nodes which will get a zero out degree to set
////		IAssert(prunable->IsNode(k));
////		TNGraph::TNodeI kiter = prunable->GetNI(k);
////		IAssert(kiter.GetOutDeg() > 0);
////		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
////			int lid = kiter.GetInNId(inEdgeNumber);
////			TNGraph::TNodeI l = prunable->GetNI(lid);
////			if (l.GetOutDeg() == 1) {
////				zeroOutDegrees[l.GetId()] = true;
////			}
////		}
////		//update the priorities of all nodes k is pointing to
////		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
////			int lid = kiter.GetOutNId(outEdgeNumber);
////			TNGraph::TNodeI l = prunable->GetNI(lid);
////			//We add one to indicate that even a node with 0 in degree is still a valid node.
////			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
////			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
////		}
////		//finalize
////		prunable->DelNode(k);
////		todo[k] = false;
////		finalOrder.Add(k);
////		if (finalOrder.Len() % infoFrequency == 0) {
////			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////		}
////
////	}
////	IAssert(startSize == finalOrder.Len());
////	IAssert(todo.find_first() == todo.npos);
////	IAssert(allNumbersIn(finalOrder));
//	return finalOrder;
//}


}

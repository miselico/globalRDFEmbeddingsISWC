/*
 * utils.cpp
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#include "utils.h"

//template<class NodeData, class EdgeData>
//TPt<TNodeEdgeNet<NodeData, EdgeData> > reverseGraph(TPt<TNodeEdgeNet<NodeData, EdgeData> > baseGraph) {
//	TPt<TNodeEdgeNet<NodeData, EdgeData> > newNet = TNodeEdgeNet<NodeData, EdgeData>::New();
//	for (typename TNodeEdgeNet<NodeData, EdgeData>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
//		newNet->AddNode(NI.GetId(), NI.GetDat());
//	}
//	//add all edges reversed
//
//	for (typename TNodeEdgeNet<NodeData, EdgeData>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
//		newNet->AddEdge(EI.GetDstNId(), EI.GetSrcNId(), EI.GetId(), EI.GetDat());
//	}
//	return newNet;
//}



TPt<TNodeEdgeNet<TStr, TStr> > reverseGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	TPt<TNodeEdgeNet<TStr, TStr> > newNet = TNodeEdgeNet<TStr, TStr>::New();
	for (typename TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}
	//add all edges reversed

	for (typename TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		newNet->AddEdge(EI.GetDstNId(), EI.GetSrcNId(), EI.GetId(), EI.GetDat());
	}
	return newNet;
}

/*
 * RDF2Walk.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "RDF2Walk.h"
#include "nTripleParser.h"

namespace {
TPt<TNodeEdgeNet<TStr, TStr> > pruneZeroOutDegree(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	//currently making a copy, no idea wheter it is smart to iterate and remove stuff at the same time
	TPt<TNodeEdgeNet<TStr, TStr> > newNet = TNodeEdgeNet<TStr, TStr>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		if (NI.GetOutDeg() > 0) {
			newNet->AddNode(NI.GetId(), NI.GetDat());
		}
	}

	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		if (newNet->IsNode(EI.GetDstNId())) {
			newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), EI.GetDat());
		}
	}
	return newNet;
}

static TStr RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static TStr OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");

}
TPair<TPt<TNodeEdgeNet<TStr, WeightedPredicate> >, TVec<TInt> > buildWalkableGraphIgnoringLiteralsAndLeafs(const TStr & filename, const GraphWeigher & weighingStrategy) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	//we need all nodes having type Thing

	TVec<TInt> entities;
	for (int i = 0; i < graph->GetNodes(); i++) {
		TNodeEdgeNet<TStr, TStr>::TNodeI candidateNode = graph->GetNI(i);
		for (int outEdgeNr = 0; outEdgeNr < candidateNode.GetOutDeg(); ++outEdgeNr) {
			TStr predicate = candidateNode.GetOutEDat(outEdgeNr);
			if (predicate == RDF_TYPE) {
				TStr object = candidateNode.GetOutNDat(outEdgeNr);
				if (object == OWL_THING) {
					entities.Add(i);
				}
			}
		}
	}

	graph = pruneZeroOutDegree(graph);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return TPair<TPt<TNodeEdgeNet<TStr, WeightedPredicate> >, TVec<TInt> >(weightedGraph, entities);
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraphIgnoringLiterals(const TStr & filename, const GraphWeigher & weighingStrategy) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return weightedGraph;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraph(const TStr & filename, const GraphWeigher & weighingStrategy) {
	TPt<TNodeEdgeNet<TStr, TStr> > graph;
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraph(filename);
	graph = graphAndNodeIndex.Val1;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return weightedGraph;
}

void TextFileSink::consume(const TVec<TStr> & path) {
	int sepCount = 0;
	for (int partNr = 0; partNr < path.Len(); partNr++) {
		fwrite(" ", sizeof(char), sepCount, fout);
		sepCount = 1;
		TStr thePart = path[partNr];
		fwrite(thePart.CStr(), sizeof(char), thePart.Len(), fout);
	}
	fwrite("\n", sizeof(char), 1, fout);
}

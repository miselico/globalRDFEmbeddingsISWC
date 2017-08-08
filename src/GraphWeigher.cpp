/*
 * GraphWeigher.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWeigher.h"

#include <iostream>

using namespace std;

namespace { //helpers

//count freq of each property:
THash<TStr, TFlt> absolute_predicate_freq(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	THash<TStr, TFlt> absolute_freq;
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr predicate = EI.GetDat();
		TFlt start = 0.0;
		absolute_freq.IsKeyGetDat(predicate, start);
		absolute_freq.AddDat(predicate, start + 1);
	}
	return absolute_freq;
}

//counts the absolute frequence of a the objects
THash<TStr, TFlt> absolute_object_freq(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	THash<TStr, TFlt> absolute_freq;
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr object = EI.GetDstNDat();
		TFlt start = 0.0;
		absolute_freq.IsKeyGetDat(object, start);
		absolute_freq.AddDat(object, start + 1);
	}
	return absolute_freq;
}

THash<TStrPr, TFlt> absolute_predicate_object_freq(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	THash<TStrPr, TFlt> absolute_freq;
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		TStr obj = EI.GetDstNDat();
		TStrPr PO(pred, obj);
		TFlt start = 0.0;
		absolute_freq.IsKeyGetDat(PO, start);
		absolute_freq.AddDat(PO, start + 1);
	}
	return absolute_freq;
}

//inverse all freq
template<class type> THash<type, TFlt> inverse_the_frequency(THash<type, TFlt> & absolute_freq) {
	THash<type, TFlt> inverse_freq;
	for (typename THash<type, TFlt>::TIter iter = absolute_freq.BegI(); iter < absolute_freq.EndI(); iter++) {
		inverse_freq.AddDat(iter.GetKey(), 1.0 / iter.GetDat());
	}
	return inverse_freq;
}

//newnetcopynodes
TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNetCopyNodes(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}
	return newNet;
}

//normalize all weights in unbalanced (sum weight on outedges == 1)
void normalize(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > unbalanced) {
	for (TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI NI = unbalanced->BegNI(); NI < unbalanced->EndNI(); NI++) {
		int node_i_outdeg = NI.GetOutDeg();
		double totalWeight = 0;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			totalWeight += NI.GetOutEDat(outEdge).W();
		}
		double totalWeightInverse = 1.0 / totalWeight;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			WeightedPredicate & wp = NI.GetOutEDat(outEdge);
			double normalized_weight = wp.W() * totalWeightInverse;
			wp.Val2 = normalized_weight;
		}
	}
}

}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > UniformWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		WeightedPredicate wpred(pred, 1);
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize(newNet);
	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > PredicateFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	THash<TStr, TFlt> absolute_freq = absolute_predicate_freq(baseGraph);
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		WeightedPredicate wpred(pred, absolute_freq.GetDat(pred));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize(newNet);
	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > InversePredicateFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	//count freq of each property:
	THash<TStr, TFlt> absolute_freq = absolute_predicate_freq(baseGraph);
	THash<TStr, TFlt> inverse_freq = inverse_the_frequency(absolute_freq);
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		WeightedPredicate wpred(pred, inverse_freq.GetDat(pred));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize(newNet);
	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > ObjectFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	THash<TStr, TFlt> absolute_freq = absolute_object_freq(baseGraph);
	PushDownWeigher subWeigher(absolute_freq, -1);
	return subWeigher.weigh(baseGraph);
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > InverseObjectFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	THash<TStr, TFlt> absolute_freq = absolute_object_freq(baseGraph);
	THash<TStr, TFlt> inverse_freq = inverse_the_frequency(absolute_freq);
	PushDownWeigher subWeigher(inverse_freq, -1);
	return subWeigher.weigh(baseGraph);
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > InverseObjectFrequencyWeigherSplitDown::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	THash<TStr, TFlt> absolute_freq = absolute_object_freq(baseGraph);
	THash<TStr, TFlt> inverse_freq = inverse_the_frequency(absolute_freq);
	SplitDownWeigher subWeigher(inverse_freq, -1);
	return subWeigher.weigh(baseGraph);
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > PredicateObjectFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	THash<TStrPr, TFlt> absolute_freq = absolute_predicate_object_freq(baseGraph);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		TStr obj = EI.GetDstNDat();
		TStrPr PO(pred, obj);
		WeightedPredicate wpred(pred, absolute_freq.GetDat(PO));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize(newNet);
	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > InversePredicateObjectFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	THash<TStrPr, TFlt> absolute_freq = absolute_predicate_object_freq(baseGraph);
	THash<TStrPr, TFlt> inverse_freq = inverse_the_frequency(absolute_freq);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr pred = EI.GetDat();
		TStr obj = EI.GetDstNDat();
		TStrPr PO(pred, obj);
		WeightedPredicate wpred(pred, inverse_freq.GetDat(PO));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize(newNet);
	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > PushDownWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);

	//add all edges with assigned weight
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		TFlt weight = this->defaultWeight;
		//overwrite default if a value is set
		nodeWeights.IsKeyGetDat(NI.GetDat(), weight);
		int node_i_indeg = NI.GetInDeg();
		if (node_i_indeg > 0) {
			if (this->defaultWeight == -1.0 && weight == -1.0) {
				cerr << "Defaultweight was  -1 and node was not found in the given map -> ERROR, quitting";
				exit(1);
			}
			//add all *in* edges with the weight
			for (int inEdge = 0; inEdge < node_i_indeg; ++inEdge) {
				TStr label = NI.GetInEDat(inEdge);
				const WeightedPredicate pred(label, weight);
				newNet->AddEdge(NI.GetInNId(inEdge), NI.GetId(), NI.GetInEId(inEdge), pred);
			}
		}
	}
	normalize(newNet);

	return newNet;

}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > SplitDownWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);

	//add all edges with assigned weight
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		TFlt weight = this->defaultWeight;
		//overwrite default if a value is set
		nodeWeights.IsKeyGetDat(NI.GetDat(), weight);
		int indeg = NI.GetInDeg();
		if (indeg > 0) {
			if (this->defaultWeight == -1.0 && weight == -1.0) {
				cerr << "Defaultweight was  -1 and node was not found in the given map -> ERROR, quitting";
				exit(1);
			}
			//add all *in* edges with the weight/indegree
			weight = weight / double(indeg);
			for (int inEdge = 0; inEdge < indeg; ++inEdge) {
				TStr label = NI.GetInEDat(inEdge);
				const WeightedPredicate pred(label, weight);
				newNet->AddEdge(NI.GetInNId(inEdge), NI.GetId(), NI.GetInEId(inEdge), pred);
			}
		}
	}
	normalize(newNet);

	return newNet;

}

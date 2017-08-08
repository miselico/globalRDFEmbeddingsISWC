/*
 * GraphWeigher.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef GRAPHWEIGHER_H_
#define GRAPHWEIGHER_H_

#include "Snap.h"
#include "WeightedPredicate.h"

class GraphWeigher {
protected:
	GraphWeigher() {

	}

public:
	virtual ~GraphWeigher() {

	}
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const = 0;
};

class UniformWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class InversePredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class PredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class ObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class InverseObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class PredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class InversePredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class InverseObjectFrequencyWeigherSplitDown: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

/**
 * Assigns to each inedge the weight assigned to the nodes.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class PushDownWeigher: public GraphWeigher {
	const THash<TStr, TFlt> nodeWeights;
	const double defaultWeight;

public:
	/**
	 * If defaultweight is set to -1, it indicates that all weights MUST be in the nodeWeights. If not, the program will be terminated.
	 */
	PushDownWeigher(const THash<TStr, TFlt> nodeWeights, const double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

/**
 * Assigns to each inedge the weight assigned to the nodes divided by the number of in edges.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight divided by #inedge
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node / #inedge
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class SplitDownWeigher: public GraphWeigher {
	const THash<TStr, TFlt> nodeWeights;
	const double defaultWeight;

public:
	/**
	 * If defaultweight is set to -1, it indicates that all weights MUST be in the nodeWeights. If not, the program will be terminated.
	 */
	SplitDownWeigher(const THash<TStr, TFlt> nodeWeights, const double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

#endif /* GRAPHWEIGHER_H_ */

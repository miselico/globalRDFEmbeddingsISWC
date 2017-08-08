/*
 * GraphWalker.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef GRAPHWALKER_H_
#define GRAPHWALKER_H_

#include "Snap.h"
#include "WeightedPredicate.h"

class GraphWalker {
protected:
	GraphWalker() {
	}

	virtual ~GraphWalker() {
	}

public:
	virtual TVec<TStr> performWalk(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph, const TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI & start) = 0;
};

/**
 * Walks on the graph using the weighted edges as probabilities for following
 *
 * Starts each walk at a randomly chosen node
 */
class RandomProportionalWalker: public GraphWalker {

private:
	const int walklength;
	TRnd random;
public:

	/**
	 * Make a walker for given length and seed to steer rng.
	 * Note, the current implementation is such that the walk consists of
	 *
	 * <startnode> <P1> node_1 <P2> node_2 ... <Pwalklength> node_walklength
	 */
	RandomProportionalWalker(int walkLength, int seed) :
			walklength(walkLength), random(TRnd(seed)) {
	}

	/**
	 * Performs a walk on the graph. Returns the walk.
	 *
	 * The walk length is at most 2*walklength+1
	 */
	virtual TVec<TStr> performWalk(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph, const TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI & startNode);

};


/*
 * Only walks of the specified length will be created. If this fail, this class will retry
 *
 * Warning: if no walks of that length can be created, this goes in an infinite loop!
 */
class LengthEnforcingWalker: public GraphWalker {
private:
	GraphWalker& actualWalker;
	const int enforcedLength;

public:
	LengthEnforcingWalker(GraphWalker& actualWalker, int enforcedLength) :
			actualWalker(actualWalker), enforcedLength(enforcedLength) {
	}
	virtual TVec<TStr> performWalk(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph, const TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI & start) {
		TVec<TStr> path;
		do {
			path = this->actualWalker.performWalk(graph, start);
		} while (path.Len() != this->enforcedLength);
		return path;
	}

};

#endif /* GRAPHWALKER_H_ */

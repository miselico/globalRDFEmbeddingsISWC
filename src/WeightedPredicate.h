/*
 * WeightedPredicate.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef WEIGHTEDPREDICATE_H_
#define WEIGHTEDPREDICATE_H_

#include "Snap.h"
//#include <cstdlib>
//#include <cmath>
//#include <cassert>
#include <limits>

class WeightedPredicate: public TPair<TStr, TLFlt> {


public:
	//Initializes weight to -INF to play safe. More likely that mistakes are found.
	WeightedPredicate() :
			TPair(TStr(), -INFINITY) {
		static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");
	}

	WeightedPredicate(TStr predicate, double weight = 1.0) :
			TPair(predicate, weight) {
	}

	const TStr P() const {
		return this->Val1;
	}

	const long double W() const {
		return this->Val2.Val;
	}

};

#endif /* WEIGHTEDPREDICATE_H_ */

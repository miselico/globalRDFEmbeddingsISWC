/*
 * utils.h
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "Snap.h"

//template<class NodeData, class EdgeData>
//TPt<TNodeEdgeNet<NodeData, EdgeData> > reverseGraph(TPt<TNodeEdgeNet<NodeData, EdgeData> > baseGraph);

TPt<TNodeEdgeNet<TStr, TStr> > reverseGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph);


static TStr RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static TStr OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");



#endif /* UTILS_H_ */

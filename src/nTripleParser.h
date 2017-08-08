/*
 * nTripleParser.h
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#ifndef NTRIPLEPARSER_H_
#define NTRIPLEPARSER_H_

#include <iostream>

#include "Snap.h" //this is a bit overkill, but includes everything needed for glib

namespace n3parser{

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraph(TStr filename);

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraphIgnoreLiterals(TStr filename);


//Triple parsetripleLine(TStr line);

TStr resourceExpander(TStr string);

}

#endif /* NTRIPLEPARSER_H_ */

/*
 * ntripleParser.cpp
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#include "nTripleParser.h"
#include "PrintTime.h"

using namespace std;

namespace {

bool DEBUGMODE = false;

class Triple: public TTriple<TStr, TStr, TStr> {
public:
	Triple(TStr S, TStr P, TStr O) :
			TTriple(S, P, O) {
	}

	TStr S() {
		return this->Val1;
	}
	TStr P() {
		return this->Val2;
	}
	TStr O() {
		return this->Val3;
	}

};

// A blank node is at the start of this line, parses it of and returns the remaining of the line.
TPair<TStr, TStr> parseBlankNodeOf(TStr line) {
	int firstUnderScoreColon = line.SearchStr("_:", 0);
	if (firstUnderScoreColon == -1) {
		cerr << "_: expected not found " << line.CStr() << endl;
	}
	TStr atBNStart = line.GetSubStr(firstUnderScoreColon);
	//We assume it is ended by space.

	TStr blankNode = line.LeftOf(' ');
	TStr rest = line.RightOf(' ');
	return TPair<TStr, TStr>(blankNode, rest);
}

TPair<TStr, TStr> parseSubjectOf(TStr line) {
	int firstUnderScoreColon = line.SearchStr("_:", 0);
	int firstAngular = line.SearchCh('<', 0);

	if (firstUnderScoreColon != -1 && firstUnderScoreColon < firstAngular) {
		return parseBlankNodeOf(line);
	}
	//otherwise it is a resource
	line = line.RightOf('<');
	TStr S = "<" + line.LeftOf('>') + ">";
	return TPair<TStr, TStr>(S, line.RightOf('>'));
}

TStr parseObject(TStr line) {

	//the line now either contains one more resource or a literal
	//Note, there could be a '<' in the string and according to production of N3 also a " in resource

	while (line[0] == ' ') {
		line = line.GetSubStr(1);
	}

	if (line.GetCh(0) == '_') {
		return parseBlankNodeOf(line).Val1;
	} else if (line.GetCh(0) == '<') {
		line = line.RightOf('<');
		return "<" + line.LeftOf('>') + ">";
	} else if (line.GetCh(0) == '"') {
		//literal
		if (DEBUGMODE && ((line.SearchStr("^^", 0) != -1) || (line.SearchStr("@", 0) != -1))) {
			cerr << "TODO language tags and datatypes are not supported properly, assumed as part of the literal : " << line.CStr() << endl;
		}
		//cut of final dot and spaces
		while (line.GetCh(line.Len() - 1) == ' ' || line.GetCh(line.Len() - 1) == '.') {
			line = line.GetSubStr(0, line.Len() - 2);
		}
		return line;
	} else {
		cerr << " Parsing error, invalid object " << line.CStr() << endl;
		exit(1);
	}

}

Triple parsetripleLine(TStr line) {
	//line = <http://www.wikidata.org/ontology#gcPrecision> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#DatatypeProperty> .
	//line = _:node1ap4j3o3kx4 <http://example.com> <http://example.com#hasname> "Test"

	TPair<TStr, TStr> S_rest = parseSubjectOf(line);

	TStr S = S_rest.Val1;

	//predicate is easy
	line = S_rest.Val2.RightOf('<');
	TStr P = line.LeftOf('>');
	P = "<" + P + ">";
	line = line.RightOf('>');

	TStr O = parseObject(line);

	return Triple(S, P, O);
}

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraphInternal(TStr filename, bool removeLiteral) {

//The graph
	TPt<TNodeEdgeNet<TStr, TStr> > Net = TNodeEdgeNet<TStr, TStr>::New();

//temporary keep track of all the nodes
	THash<TStr, int> addedNodes;

	PSIn FInPt = TFIn::New(filename);
	TStr line;
	int count = 0;
//To reuse the strings in the node/edge labels
	THashSet<TStr> stringpool;

	while (FInPt->GetNextLn(line)) {
		if (line.IsWs()) {
			continue;
		}
		if (line.SearchCh('#', 0) == 0) {
			//comment
			continue;
		}
		count++;
		if (count % 1000000 == 0) {
			cout << currentTime() << "Read " << count << " triples" << endl;
		}
		Triple valuesnotPooled = parsetripleLine(line);
		//This saves about 10% memory on a small test. Perhaps more on a larger one.
		Triple values(valuesnotPooled.S(), stringpool.GetKey(stringpool.AddKey(valuesnotPooled.P())), valuesnotPooled.O());
		if (removeLiteral && values.O().GetCh(0) == '"') {
			//do not add this literal
			continue;
		}

		int subjectIndex = 0;
		bool hasEdgePossibility = true;
		if (addedNodes.IsKeyGetDat(values.S(), subjectIndex)) {
			//nothing, subjectIndex now contains the node ID
		} else {
			//add new Node, save index
			subjectIndex = Net->AddNode(-1, values.S());
			addedNodes.AddDat(values.S(), subjectIndex);
			hasEdgePossibility = false;
		}

		int objectIndex = 0;
		if (addedNodes.IsKeyGetDat(values.O(), objectIndex)) {
			//nothing, objectIndex now contains the node ID
		} else {
			//add new Node, save index
			objectIndex = Net->AddNode(-1, values.O());
			addedNodes.AddDat(values.O(), objectIndex);
			hasEdgePossibility = false;
		}
		bool hasEdge = false;
		if (hasEdgePossibility) {
			TNodeEdgeNet<TStr, TStr>::TNodeI ni = Net->GetNI(subjectIndex);
			for (int i = 0; i < ni.GetOutDeg(); ++i) {
				int otherObjectindex = ni.GetOutNId(i);
				if (otherObjectindex == objectIndex) {
					TStr pred = ni.GetOutEDat(i);
					if (pred == values.P()) {
						hasEdge = true;
						break;
					}
				}
			}
		}
		if (!hasEdge) {

			//add edge
			Net->AddEdge(subjectIndex, objectIndex, -1, values.P());
		}

	}

	stringpool.Clr(true, -1);

	addedNodes.Pack();

	return TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> >(Net, addedNodes);
}

}

namespace n3parser {

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraph(TStr filename) {
	return buildRDFGraphInternal(filename, false);
}

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraphIgnoreLiterals(TStr filename) {
	return buildRDFGraphInternal(filename, true);

}
}

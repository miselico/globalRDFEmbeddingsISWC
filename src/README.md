

This repository contains code related to embedding entities in RDF graphs into a vector space.


Overview of files:
=====================


BCA
----
Bookmark coloring algorithm, and pushed-Bookmark coloring algorithm.
Implementation of approximate personal pagerank

doublePriorityQueue
--------------------
The priorityQueueu in glib uses floats, (weirdly converting to doubles internally).
This provides a priority queue using doubles for the priority.
This is used in BCA

GraphWalker
-------------
Strategies for walking over a graph where edges have a weight

GraphWeigher
-------------
Converters adding weights to a directed, labeled graph.
Other conversions could be done by these 'weighers' as well. The final naming for these might need some thinking.

Main
------
The main() has to be somewhere where it can be found...

nTriplesParser
---------------
A quick 'n dirty parser for n-triples. No error handling, very sensitive to errors in input.

RandomWalkExperiments
---------------------
Experiments the generation of random walks

RDF2Co_occurence
-------------------
The idea is to create co-occurence matrices from RDF graphs. This is not yet completely integrated with glove.

RDF2Walk
-----------
Helper functions for going from a n-triples file to a file containing walks

WeigtedPredicate
----------------
The type of the value used on the arc of a weighted graph. Includes the weight and a predicate




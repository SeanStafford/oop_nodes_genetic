#pragma once
#include <bitset>
#include <vector>
#include <map>
#include "CompileTimeConstants.h"

using namespace std;


class Node {
public:
	Node(int& former_highest_id) :id(++former_highest_id) {
		recently_updated = true;
	}

	void PickEdgesToForm(list<Node>&, list<Node>::iterator, vector<double>, vector<double>, vector<bool>, vector<bool>);

	void FormInwardEdge(Node*, float);

	void FormOutwardEdge(Node*, float);

	void UpdateSum();

	void Mutate(bitset<genome_size>, uniform_real_distribution<double>, default_random_engine);

	double ReturnSum() { return fitness; }

	bool CheckRecentlyUpdated() { return recently_updated; }

	void ResetRecentlyUpdated() { recently_updated = false; }

	void Disentangle();

	int CountInDegree() { return inwardly_directed_edges.size(); }

	int CountOutDegree() { return outwardly_directed_edges.size(); }

	int ReturnID() { return id; }

	vector<double> ReturnOutWardEdges(list<Node>&, list<Node>::iterator);

	bitset<genome_size> ReturnGenome() { return genome; }

	void AssignGenome(int gen) { genome = gen; }

protected:
	// Flag for whether UpdateNodeList can skip this node
	bool recently_updated;
	
	// Counts how many nodes have been created as of the creation of this node 
	int id;

	// Sum of weights of all inwardly directed edges
	double fitness; 

	// Maps Node pointers to interaction coefficients of nodes affecting this node's fitness
	map<Node*, double> inwardly_directed_edges;
	// Maps Node pointers to interaction coefficients of nodes whose fitness is affected by this node
	map<Node*, double> outwardly_directed_edges;

	// Bitstring that represents this node's genome
	bitset<genome_size> genome;
};

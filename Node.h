#pragma once
#include <bitset>
#include <vector>
#include <map>
#include <list>
#include <random>

using namespace std;

typedef bitset<64> genome_t;

class Node {
public:
	static int highest_id;
	Node() :id(++highest_id) {
		recently_updated = true;
	}

	void PickEdgesToForm(list<Node>&, const vector<double>& X, const vector<double>& Y, const vector<bool>& cX, const vector<bool>& cY);

	void FormInwardEdge(Node*, float);

	void FormOutwardEdge(Node*, float);

	void UpdateSum();

	void Mutate(genome_t, uniform_real_distribution<double>, default_random_engine, long genome_size);

	double ReturnSum() { return fitness; }

	bool CheckRecentlyUpdated() { return recently_updated; }

	void ResetRecentlyUpdated() { recently_updated = false; }

	void Disentangle();

	int CountInDegree() { return inwardly_directed_edges.size(); }

	int CountOutDegree() { return outwardly_directed_edges.size(); }

	int ReturnID() const { return id; }

	genome_t ReturnGenome() const { return genome; }

	void AssignGenome(int gen) { genome = gen; }

protected:
	// Flag for whether UpdateNodeList can skip this node
	bool recently_updated;
	
	// Counts how many nodes have been created as of the creation of this node 
	const int id;

	// Sum of weights of all inwardly directed edges
	double fitness; 

	// Maps Node pointers to interaction coefficients of nodes affecting this node's fitness
	map<Node*, double> inwardly_directed_edges;
	// Maps Node pointers to interaction coefficients of nodes whose fitness is affected by this node
	map<Node*, double> outwardly_directed_edges;

	// Bitstring that represents this node's genome
	genome_t genome;
};

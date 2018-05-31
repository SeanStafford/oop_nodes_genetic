#include "Node.h"
#include "Functions.h"

int Node::highest_id = 0;
int Node::death_count = 0;
int Node::burn_in_id;
vector<int> Node::lifespans;


// Randomly choose one bit to flip. Flip the choosen bit. Then set the modified bitstring as the genome.
void Node::Mutate(genome_t old_genome, uniform_real_distribution<double> dist, default_random_engine gen, long genome_size) {
	int rand_index = dist(gen) * genome_size;
	old_genome[rand_index] = !old_genome[rand_index];
	genome = old_genome;
}

// Forms edges between a newly created node and existing nodes
void Node::PickEdgesToForm(list<Node>& node_list, const vector<double>& X, const vector<double>& Y, const vector<bool>& connect_x, const vector<bool>& connect_y) {
	// Cycle through all nodes except this one
	for( Node& node : node_list ) {
		if(&node == this) {
			continue;
		}

		//// Find the interaction coefficient from the current node to the node i + 1 along in the list. If the coefficient
		//// is nonzero form store this information in the respective pointer maps of the relevant nodes
		double coeff = GetCoefficient(genome, node.ReturnGenome(), X, Y, connect_x, connect_y);
		if (coeff) { FormOutwardEdge(&node, coeff); }

		//// Find the interaction coefficient from the node i + 1 along in the list to the current node. If the coefficient
		//// is nonzero form store this information in the respective pointer maps of the relevant nodes
		coeff = GetCoefficient(node.ReturnGenome(), genome, X, Y, connect_x, connect_y);
		if (coeff) { FormInwardEdge(&node, coeff); }

		//// Note: While the initial community is being formed, running both FormOutwardEdge and FormInwardEdge is redundant.
		//// However, this is fine, because the map insertion operation does nothing if you attempt to dublicate keys.
	}
}

void Node::FormOutwardEdge(Node* node_receiving_edge, double coeff) {
	// Store a pointer to the node whose fitness you are affecting along with the coefficient of interaction in your outward map
	outwardly_directed_edges.insert(pair<Node*, double>(node_receiving_edge, coeff));
	// Go to the the node whose fitness you are affecting. Store a pointer to yourself along with the coefficient of interaction in
	// the node's inward map.
	(node_receiving_edge->inwardly_directed_edges).insert(pair<Node*, double>(this, coeff));
	// If the new edge has a negative weight, it may make fitness nonpositive so the affected node is flagged
	if (coeff <= 0) { node_receiving_edge->recently_updated = true; }
}

void Node::FormInwardEdge(Node* node_giving_edge, double coeff) {
	// Store a pointer to the node which is affecting your fitness along with the coefficient of interaction in your inward map
	inwardly_directed_edges.insert(pair<Node*, double>(node_giving_edge, coeff));
	// If the new edge has a negative weight, it may make fitness nonpositive so this node is flagged
	if (coeff <= 0) { recently_updated = true; }
	// Go to the node which is affecting your fitness. Store a pointer to yourself along with the coefficient of interaction in
	// the node's outward map.
	(node_giving_edge->outwardly_directed_edges).insert(pair<Node*, double>(this, coeff));
}

void Node::UpdateSum() {
  double sum = 0.0;
	for( auto& pair : inwardly_directed_edges ) {
		sum += pair.second;
	}
	fitness = sum;
}

// If I do not remove the map entries of species that have died, the program crashes. This method takes care of that
void Node::Disentangle() {
	for(auto & inlink : inwardly_directed_edges) {
		Node* neighbor = inlink.first;
		auto& n_links = neighbor->outwardly_directed_edges;
		n_links.erase( n_links.find(this) );
	}
  for(auto& outlink: outwardly_directed_edges) {
    Node* neighbor = outlink.first;
		if(outlink.second > 0.0) { neighbor->recently_updated = true; }
		auto& n_links = neighbor->inwardly_directed_edges;
		n_links.erase( n_links.find(this) );
	}
}

// Keep track of all lifesspans in order to plot power spectrum of lifetimes later
void Node::InitiateLifespanDistribution(int initial_size, int total_steps, int burn_in) {
	lifespans.resize(initial_size + total_steps);
	burn_in_id = burn_in + initial_size;
}

void Node::UpdateLifespanDistribution() {
	cout << "Reached UpdateLifespanDistribution" << endl;
	cout << "highest_id = " << highest_id;
	cout << " burn_in_id = " << burn_in_id << endl;
	if (highest_id >= burn_in_id) {
		death_count++;
		lifespans[death_count] = highest_id - id;
	}
	cout << "Left UpdateLifespanDistribution" << endl;
}
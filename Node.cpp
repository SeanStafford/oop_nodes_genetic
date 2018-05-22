#include "Node.h"
#include "Functions.h"

int Node::highest_id = 0;

// Randomly choose one bit to flip. Flip the choosen bit. Then set the modified bitstring as the genome.
void Node::Mutate(genome_t old_genome, uniform_real_distribution<double> dist, default_random_engine gen, long genome_size) {
	int rand_index = dist(gen) * genome_size;
	old_genome[rand_index] = !old_genome[rand_index];
	genome = old_genome;
}

// Forms edges between a newly created node and existing nodes
void Node::PickEdgesToForm(list<Node>& node_list, list<Node>::iterator this_element, vector<double> X, vector<double> Y, vector<bool> connect_x, vector<bool> connect_y) {
	// Cycle through all nodes except this one
	for (int i = 0; i < node_list.size() - 1; ++i) {
		IterateCyclically(this_element, node_list);
		//// Find the interaction coefficient from the current node to the node i + 1 along in the list. If the coefficient
		//// is nonzero form store this information in the respective pointer maps of the relevant nodes
		float coeff = GetCoefficient(genome, this_element->ReturnGenome(), X, Y, connect_x, connect_y);
		if (coeff) { FormOutwardEdge(&(*this_element), coeff); }

		//// Find the interaction coefficient from the node i + 1 along in the list to the current node. If the coefficient
		//// is nonzero form store this information in the respective pointer maps of the relevant nodes
		coeff = GetCoefficient(this_element->ReturnGenome(), genome, X, Y, connect_x, connect_y);
		if (coeff) { FormInwardEdge(&(*this_element), coeff); }

		//// Note: While the initial community is being formed, running both FormOutwardEdge and FormInwardEdge is redundant.
		//// However, this is fine, because the map insertion operation does nothing if you attempt to dublicate keys.
	}
}

void Node::FormOutwardEdge(Node* node_receiving_edge, float coeff) {
	// Store a pointer to the node whose fitness you are affecting along with the coefficient of interaction in your outward map
	outwardly_directed_edges.insert(pair<Node*, double>(node_receiving_edge, coeff));
	// Go to the the node whose fitness you are affecting. Store a pointer to yourself along with the coefficient of interaction in
	// the node's inward map.
	(node_receiving_edge->inwardly_directed_edges).insert(pair<Node*, double>(this, coeff));
	// If the new edge has a negative weight, it may make fitness nonpositive so the affected node is flagged
	if (coeff <= 0) { node_receiving_edge->recently_updated = true; }
}

void Node::FormInwardEdge(Node* node_giving_edge, float coeff) {
	// Store a pointer to the node which is affecting your fitness along with the coefficient of interaction in your inward map
	inwardly_directed_edges.insert(pair<Node*, double>(node_giving_edge, coeff));
	// If the new edge has a negative weight, it may make fitness nonpositive so this node is flagged
	if (coeff <= 0) { recently_updated = true; }
	// Go to the node which is affecting your fitness. Store a pointer to yourself along with the coefficient of interaction in
	// the node's outward map.
	(node_giving_edge->outwardly_directed_edges).insert(pair<Node*, double>(this, coeff));
}

void Node::UpdateSum() {
	fitness = 0;
	map<Node*, double>::iterator temp_itr;
	for (temp_itr = inwardly_directed_edges.begin(); temp_itr != inwardly_directed_edges.end(); ++temp_itr) {
		fitness += temp_itr->second;
	}
}

// If I do not remove the map entries of species that have died, the program crashes. This method takes care of that
void Node::Disentangle() {
	map<Node*, double>::iterator temp_itr_in;
	for (temp_itr_in = inwardly_directed_edges.begin(); temp_itr_in != inwardly_directed_edges.end(); ++temp_itr_in) {
		map<Node*, double>& out_ptrs = (temp_itr_in->first)->outwardly_directed_edges;
		out_ptrs.erase(out_ptrs.find(this));
	}
	map<Node*, double>::iterator temp_itr_out;
	for (temp_itr_out = outwardly_directed_edges.begin(); temp_itr_out != outwardly_directed_edges.end(); ++temp_itr_out) {
		// Removing a positive edge to a node could cause its fitness to become non-positive so I flag it when this happens
		if (temp_itr_out->second > 0) { (temp_itr_out->first)->recently_updated = true; }
		map<Node*, double>& in_ptrs = (temp_itr_out->first)->inwardly_directed_edges;
		in_ptrs.erase(in_ptrs.find(this));
	}
}

// This method is left over from when I was using the migration model and I have not used it recently
vector<double> Node::ReturnOutWardEdges(list<Node>& node_list, list<Node>::iterator this_element) {

	int node_tally = node_list.size();

	int initial_size = outwardly_directed_edges.size();

	vector<double> alloutwardedges(node_tally, 0);

	vector<vector<double>> outwardedgeinfo(initial_size, vector<double>(2));

	int i = 0;
	map<Node*, double>::iterator temp_itr_out;
	for (temp_itr_out = outwardly_directed_edges.begin();
		temp_itr_out != outwardly_directed_edges.end(); ++temp_itr_out) {
		outwardedgeinfo[i][0] = (temp_itr_out->first)->id;
		outwardedgeinfo[i][1] = temp_itr_out->second;
		i++;
	}

	list<Node>::iterator itr = node_list.begin();

	int next_outward_edge_index = 0;
	for (int i = 1; i < node_tally; ++i) {
		IterateCyclically(itr, node_list);
		if (next_outward_edge_index >= outwardedgeinfo.size()) { break; }
		else if (itr->id == outwardedgeinfo[next_outward_edge_index][0]) {
			alloutwardedges[i] = outwardedgeinfo[next_outward_edge_index][1];
			next_outward_edge_index++;
		}
	}
	return alloutwardedges;
}

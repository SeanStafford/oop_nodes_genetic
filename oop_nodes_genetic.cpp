//#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <list>
#include <vector>
#include <map>
#include <fstream>
#include <random>
#include <string>
#include <time.h>
#include <bitset>
#include <math.h>
#include "Node.h"
#include "Functions.h"
#include "CompileTimeConstants.h"



using namespace std;
int main(int argc, char* argv[]) {
	// First I put all command line parameters into an array
	//// I declare the array
	double model_parameters[num_input_params];
	//// I modify the array within the function to include command line parameters
	GetModelParameters(model_parameters, argc, argv);

	// Initialize random generators
	// Use seed given from command line unless given seed is 0. Then uses time(0) to generate seed randomly
	//// I initialize default_random_generator and a normal_distribution and a uniform_real_distribution
	//// I use the normal_distribution to generate the interaction coefficients
	//// I use the uniform_real_distribution for the following:
	////// generate the interaction bools
	////// choose which species are in the initial community
	////// choose which species is mutated to form new species
	////// choose which bit is flipped when creating a mutant species
	default_random_engine generator(model_parameters[8] ? model_parameters[8] : time(0));
	normal_distribution<double> norm_dist(0, 1);
	uniform_real_distribution<double> unif_dist(0.0, 1.0);

	// I initialize the X, Y, connect_x and connect_y vectors that will be used to generate the interaction matrix
	vector<double> X(3 * genome_space, 0);
	vector<double> Y(3 * genome_space, 0);
	vector<bool> connect_x(3 * genome_space, false);
	vector<bool> connect_y(3 * genome_space, false);
	//// I populate X & Y using the default_random_generator and normal_distribution
	for (int i = 0; i < X.size(); i++) { X[i] = norm_dist(generator); }
	for (int i = 0; i < Y.size(); i++) { Y[i] = norm_dist(generator); }
	//// I populate connect_x and connect_y using uniform_real_distribution to determine if the bool should be true
	//// The probability of connection, c, is stored as model_parameters[5]. The AND operator is used with a bool from each 
	//// vector to determine if they are connected to the probability of a bool in each vector being true must be
	//// sqrt(c) to make the probability of connection c
	for (int i = 0; i < connect_x.size(); i++) {
		if (unif_dist(generator) < sqrt(model_parameters[5])) { connect_x[i] = true; }
	}
	for (int i = 0; i< connect_y.size(); i++) {
		if (unif_dist(generator) < sqrt(model_parameters[5])) { connect_y[i] = true; }
	}

	// Sometimes I am just interested in checking how correlated the interaction matrix is.
	// If so I will have specified so in the commandline parameters.
	// I just run CheckCorrelations and end the program
	if (model_parameters[0]) {
		CheckCorrelations(X, Y, connect_x, connect_y);
		return 0;
	}


	// Initialize linked list for Node objects
	list<Node> node_list;
	//// A counter is used to count how many total nodes have been made (for debugging purposes)
	//// The counter is advanced by 1 in the Node Constructor function
	int highest_id = 0;
	//// The linked list is filled with Node objects that don't have any attributes yet except the node counter
	//// model_parameters[4] is the commandline parameter for initial community size 
	for (int i = 0; i < model_parameters[4]; ++i) { node_list.push_back(Node(highest_id)); }
	//// cout a progress report
	PrintOutMessage(1);

	// Randomly assign unique genomes to each Node object
	//// GenerateInitialGenomes creates a vector of the decimal representations of choosen genomes
	vector<int> initial_genomes = GenerateInitialGenomes(model_parameters[4], unif_dist, generator);
	//// Initialize iterator for vector
	list<Node>::iterator itr = node_list.begin();
	for (int i = 0; i < node_list.size(); ++i) {
		itr->AssignGenome(initial_genomes[i]);
		IterateCyclically(itr, node_list);
	}
	//// cout a progress report
	PrintOutMessage(2);

	// Now iterate through intitial community, generating edges for each new node
	for (int i = 0; i < node_list.size(); ++i) {
		itr->PickEdgesToForm(node_list, itr, X, Y, connect_x, connect_y);
		IterateCyclically(itr, node_list);
	}
	//// cout a progress report
	PrintOutMessage(3);

	// Dispose of any Node objects whose species do not survive
	//// model_parameters[2] is a bool from commandline that signals whether 0 fitness species survive
	//// UpdateNodeList only returns False, ending the while loop, once all remaining species survive
	while (UpdateNodeList(node_list, model_parameters[2])) {};
	//// cout a progress report
	PrintOutMessage(4);
	cout << "The initial population size is " << node_list.size() << endl;

	// Create output file for time series data
	ofstream output_file("TimeSeries.txt");

	// Main While Loop
	//// Initialize a timestep counter for the main while loop
	int step = 0;
	//// While loop can be ended if either (1) all Node objects are gone or (2) The number of steps has exceeded the
	//// the amount specified by the commandline parameter model_parameters[6]
	while (node_list.size() && (++step <= model_parameters[6])) {
		// Add a new mutant species to the network
		//// First, use uniform_real_distribution to select one of the existing species. The mutant introduced to the network will be
		//// mutated from this species
		int temp_rand_index = unif_dist(generator) * node_list.size();
		//// Reset iterator because it may be pointing to a Node object that was eliminated in an UpdateNodeList step
		//// Then set iterator to selected species and retrieve its genome
		itr = node_list.begin();
		IterateCyclically(itr, node_list, temp_rand_index);
		bitset<genome_size> old_genome = itr->ReturnGenome();
		//// Now add a species without specifying its genome and set the iterator to the newly added species.
		node_list.push_back(Node(highest_id));
		itr = --node_list.end();
		//// Then use Mutate to set the species' genome to a bitstring 1 Hamming distance from old_genome
		itr->Mutate(old_genome, unif_dist, generator);

		// Check if mutant already existed among living species
		//// Retrieve the newly created genome
		bitset<genome_size> new_genome = itr->ReturnGenome();
		//// Create a flag for detecting whether the genome is a duplicate
		bool new_is_duplicate = false;
		//// Check each nodes to see if its genome matches new_genome until (1) the duplicate flag is set off or (2) all
		//// old nodes have been checked (the new one will always match itself obviously)
		for (int i = 0; i < node_list.size() - 1; ++i) {
			if (!new_is_duplicate) {
				IterateCyclically(itr, node_list);
				if (itr->ReturnGenome() == new_genome) { new_is_duplicate = true; }
			}
		}

		// Update node according to whether the genome was a duplicate
		//// Reset iterator to new node
		itr = --node_list.end();
		//// If it is a duplicate delete it
		if (new_is_duplicate) { node_list.erase(itr); }
		//// Otherwise generate edges for each node and dispose of any nodes that can no longer survive
		else {
			itr->PickEdgesToForm(node_list, itr, X, Y, connect_x, connect_y);
			while (UpdateNodeList(node_list, model_parameters[2])) {};
		}

		// Every x steps print the network size to the output file. x is specified by the commandline parameter model_parameters[7]
		if (step % (int)model_parameters[7] == 0) {
			output_file << node_list.size() << endl;
		}
	}

	// cout the final population so I check the stdout file to see if the network survived.
	cout << "The final population size is " << node_list.size() << endl;

	// Close output file and end program
	output_file.close();
	return 0;
}
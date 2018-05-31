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
#include "ModelParameters.hpp"
#include "Results.hpp"



using namespace std;
int main(int argc, char* argv[]) {
	// First I put all command line parameters into an array
	//// I declare the array
	ModelParameters param = ModelParameters::Load(argc, argv);
	param.Print();
	//// I modify the array within the function to include command line parameters

	// Initialize random generators
	// Use seed given from command line unless given seed is 0. Then uses time(0) to generate seed randomly
	//// I initialize default_random_generator and a normal_distribution and a uniform_real_distribution
	//// I use the normal_distribution to generate the interaction coefficients
	//// I use the uniform_real_distribution for the following:
	////// generate the interaction bools
	////// choose which species are in the initial community
	////// choose which species is mutated to form new species
	////// choose which bit is flipped when creating a mutant species
	default_random_engine generator(param.seed ? param.seed : time(0));
	normal_distribution<double> norm_dist(param.mean, 1);
	uniform_real_distribution<double> unif_dist(0.0, 1.0);

	// I initialize the X, Y, connect_x and connect_y vectors that will be used to generate the interaction matrix
	vector<double> X(3 * param.GenomeSpace(), 0);
	vector<double> Y(3 * param.GenomeSpace(), 0);
	vector<bool> connect_x(3 * param.GenomeSpace(), false);
	vector<bool> connect_y(3 * param.GenomeSpace(), false);
	//// I populate X & Y using the default_random_generator and normal_distribution
	for (int i = 0; i < X.size(); i++) { X[i] = norm_dist(generator); }
	for (int i = 0; i < Y.size(); i++) { Y[i] = norm_dist(generator); }
	//// I populate connect_x and connect_y using uniform_real_distribution to determine if the bool should be true
	//// The probability of connection, c, is stored as model_parameters[5]. The AND operator is used with a bool from each 
	//// vector to determine if they are connected to the probability of a bool in each vector being true must be
	//// sqrt(c) to make the probability of connection c
	for (int i = 0; i < connect_x.size(); i++) {
		if (unif_dist(generator) < sqrt(param.link_density)) { connect_x[i] = true; }
	}
	for (int i = 0; i< connect_y.size(); i++) {
		if (unif_dist(generator) < sqrt(param.link_density)) { connect_y[i] = true; }
	}

	// Sometimes I am just interested in checking how correlated the interaction matrix is.
	// If so I will have specified so in the commandline parameters.
	// I just run CheckCorrelations and end the program
	if (param.checking_correlation) {
		CheckCorrelations(X, Y, connect_x, connect_y, param.genome_length);
		return 0;
	}


	// Initialize linked list for Node objects
	list<Node> node_list;

	//// The linked list is filled with Node objects that don't have any attributes yet except the node counter
	//// model_parameters[4] is the commandline parameter for initial community size 
	for (int i = 0; i < param.initial_population_size; ++i) { node_list.push_back(Node()); }
	//// cout a progress report
	PrintOutMessage(1);

	// Randomly assign unique genomes to each Node object
	//// GenerateInitialGenomes creates a vector of the decimal representations of choosen genomes
	vector<int> initial_genomes = GenerateInitialGenomes(param.initial_population_size, unif_dist, generator, param.GenomeSpace());
	//// Initialize iterator for vector

    auto itr = node_list.begin();
	// Initiate lifespan distribution with size and burn in information
	int burn_in = 1000;
	itr->InitiateLifespanDistribution(param.initial_population_size, param.total_steps, burn_in);
    for (int i = 0; i < node_list.size(); ++i, ++itr) {
      itr->AssignGenome(initial_genomes[i]);
    }
	
	//// cout a progress report
	PrintOutMessage(2);

	// Now iterate through intitial community, generating edges for each new node
  for(Node& node : node_list) {
    node.PickEdgesToForm(node_list, X, Y, connect_x, connect_y);
  }
	//// cout a progress report
	PrintOutMessage(3);

	// Dispose of any Node objects whose species do not survive
	//// model_parameters[2] is a bool from commandline that signals whether 0 fitness species survive
	//// UpdateNodeList only returns False, ending the while loop, once all remaining species survive
	while (UpdateNodeList(node_list, param.zero_fitness_extinction)) {};
	//// cout a progress report
	PrintOutMessage(4);

  Results results;
  results.initial_population_size = node_list.size();
	cerr << "The initial population size is " << node_list.size() << endl;

	// Create output file for time series data and lifespan distribution data
	ofstream timeseries_file("TimeSeries.txt");
	ofstream lifespan_fle("Lifespans.txt");

	// Main While Loop
	//// Initialize a timestep counter for the main while loop
	//// While loop can be ended if either (1) all Node objects are gone or (2) The number of steps has exceeded the
	//// the amount specified by the commandline parameter model_parameters[6]
	int how_many_dublicates = 0;
  long average = 0;
  long total_steps = 0;
	for( int step = 0; step < param.total_steps; step++) {
		if( node_list.size() == 0 ) { break; }
		// Add a new mutant species to the network
		//// First, use uniform_real_distribution to select one of the existing species. The mutant introduced to the network will be
		//// mutated from this species
		int temp_rand_index = unif_dist(generator) * node_list.size();
		//// Reset iterator because it may be pointing to a Node object that was eliminated in an UpdateNodeList step
		//// Then set iterator to selected species and retrieve its genome
		auto itr = node_list.begin();
		for(int i=0; i < temp_rand_index; i++) { ++itr; }
		genome_t old_genome = itr->ReturnGenome();
		//// Now add a species without specifying its genome and set the iterator to the newly added species.
		node_list.push_back(Node());
    const auto new_node = --node_list.end();
		//// Then use Mutate to set the species' genome to a bitstring 1 Hamming distance from old_genome
		new_node->Mutate(old_genome, unif_dist, generator, param.genome_length);

		// Check if mutant already existed among living species
		//// Retrieve the newly created genome
		genome_t new_genome = new_node->ReturnGenome();
		//// Create a flag for detecting whether the genome is a duplicate
		//// Check each nodes to see if its genome matches new_genome until (1) the duplicate flag is set off or (2) all
		//// old nodes have been checked (the new one will always match itself obviously)
		bool new_is_duplicate = false;
		for( const Node& n : node_list ) {
			if( &n == &(*new_node) ) { continue; }
			if( n.ReturnGenome() == new_genome ) {
				new_is_duplicate = true;
				break;
			}
		}

		// Update node according to whether the genome was a duplicate
		//// If it is a duplicate delete it

		if (new_is_duplicate) {
			node_list.erase(new_node);
			cout << "There have been " << how_many_dublicates++ << " duplicates" << endl;
		}
		//// Otherwise generate edges for each node and dispose of any nodes that can no longer survive
		else {
			new_node->PickEdgesToForm(node_list, X, Y, connect_x, connect_y);
			while (UpdateNodeList(node_list, param.zero_fitness_extinction)) {};
		}

		// Every x steps print the network size to the output file. x is specified by the commandline parameter model_parameters[7]
		if (step % param.output_steps == 0) {
      double ave_deg = 0.0;
      for(const Node& n: node_list) { ave_deg += n.CountInDegree(); }
			timeseries_file << step << ' ' << node_list.size() << ' ' << ave_deg/node_list.size() << endl;
		}
    total_steps = step+1;
    average += node_list.size();
	}
  results.average_population_size = (double)average / total_steps;
  results.total_steps = total_steps;
  results.final_population_size = node_list.size();

	// cout the final population so I check the stdout file to see if the network survived.
	cerr << "The final population size is " << results.final_population_size << endl;

  std::ofstream fout("_output.json");
  results.PrintJson(fout);
  fout.close();

	timeseries_file.close();

	auto iter = node_list.begin();
	iter->lifespans.resize(iter->death_count);
	for (int lifespan : iter->lifespans) {
		lifespan_fle << lifespan << endl;
	}
	lifespan_fle.close();

	return 0;
}
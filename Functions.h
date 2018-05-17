#pragma once

#include <random>
#include "FunctionDeclarations.h"
#include "NodeMethods.h"

class Node;


void PrintOutMessage(int message_id) {
	if (message_id == 1) { cout << "generated nodes of initial community" << endl; }
	else if (message_id == 2) { cout << "assigned genomes to initial community" << endl; }
	else if (message_id == 3) { cout << "added edges to initial community" << endl; }
	else if (message_id == 4) { cout << "updated initial community" << endl; }
}

// This function allows me to include all the commandline parameters in the name of a file. Since I have been using OACIS I do not need it
string AssignFileName(double(&model_parameters)[num_input_params]) {
	ostringstream sstream("oop_nodes_output", ios_base::app | ios_base::out);
	for (auto param : model_parameters) { sstream << "_" << param; }
	sstream << ".txt";
	return sstream.str();
}

// This function asks for commandline parameters into an array, makes sure the input is the right variable type,
// then places them into an array. If they are the wrong type, it stops the program and asks for the right type.
void GetModelParameters(double(&model_parameters)[num_input_params], int num_cmdline_inputs, char** cmdline_input) {
	string instructions[num_input_params] = {
		"Enter 1 if checking correlations or 0 for normal run: ",
		"Enter 1 for fixed degree or 0 otherwise: ",
		"Enter 1 if zero fitness sufficient for death or 0 otherwise: ",
		"Enter 1 if edge direction is determined by a coin flip: ",
		"Enter initial population size (at most " + to_string(genome_space) + " ): ",
		"Enter fixed degree or degree density: ",
		"Enter total number of steps to perform: ",
		"Enter output step size: ",
		"Enter 0 for randomly generated seed or other integer to choose one: "
	};

	string received_cmdline = " (Received input from commandline argument.)\n";
	string input_problem = "Problem with input. Try again: ";

	string input_str;

	bool a_bool;
	int an_int;
	double a_double;

	string var_type = "bool";

	for (int i = 1; i <= num_input_params; i++) {
		cout << instructions[i - 1];
		bool first_attempt = true;

		// In order for the program to function properly, input parameters must be interpreted
		// as the follwing variable types 
		if (i <= 4) { var_type = "bool"; }
		else if (i == 5 || i == 7 || i == 8 || i == 9 || (i == 6 && model_parameters[1])) { var_type = "int"; }
		else if (i == 6) { var_type = "double"; }

		do {
			// If there is an available cmdline input and I have not already tried this, assign it to input_str.
			// Else, ask for an input to which I will assign input_str.
			if (first_attempt && (num_cmdline_inputs > i)) {
				// cmdline_input[0] is the name of the program
				input_str = cmdline_input[i];
				cout << input_str << received_cmdline;
				first_attempt = false;
			}
			else { getline(cin, input_str); }
		} while (
			// If I'm looking for an X-type variable, attempt to cast input_str to the variable I have declared
			// of X-type.
			/// If there is no casting error, the X-type variable previously declared now has the input_str as its
			/// value. The argument to the left of && evaluates false and so there is no need to evaluate the second
			/// argument and input_problem is never printed
			///
			/// If this casting experiences an error, the argument to the left of && evaluates true and so
			/// it evaluates the second argument, printing input_problem and evaluating true. This means the whole
			/// while loop evaluates to true and the do part of the loop must be repeated.
			((i == 5) ? (!(istringstream{ input_str } >> an_int) || (an_int > genome_space))
				: (var_type == "bool") ? !(istringstream{ input_str } >> a_bool)
				: ((var_type == "int") ? !(istringstream{ input_str } >> an_int)
					: !(istringstream{ input_str } >> a_double))) &&
					(cout << input_problem)
			);

		// Assign whichever variable is of the type I am looking for to the next spot in model_parameters.
		model_parameters[i - 1] = ((var_type == "bool") ? a_bool : ((var_type == "int") ? an_int : a_double));
	}
}

// This function is left over from when I was using the migration model and I have not used it recently
vector<vector<float>> MeasureDegreeDistribution(list<Node>& node_list) {
	// I will create a list with elements f_i = percent of nodes with i degree
	// Initially I will make it size 37, but it can expand if necessary
	int initial_size = 37;

	vector<vector<float>> degree_distributions(2, vector<float>(initial_size));
	int degree_value[2];

	list<Node>::iterator itr;
	for (itr = node_list.begin(); itr != node_list.end(); itr++) {
		degree_value[0] = itr->CountInDegree();
		degree_value[1] = itr->CountOutDegree();

		for (int i = 0; i <= 1; i++) {
			if (degree_value[i] > degree_distributions[i].size()) {
				degree_distributions[i].resize(degree_value[i] + 1);
			}
			degree_distributions[i][degree_value[i]]++;
		}
	}

	int node_tally = node_list.size();

	for (int i = 0; i <= 1; i++) {
		for (int j = 0; j < degree_distributions[i].size(); j++) {
			degree_distributions[i][j] = degree_distributions[i][j] / node_tally;
		}
	}

	return degree_distributions;
}

// This function is left over from when I was using the migration model and I have not used it recently
vector<vector<double>> RetrieveGraphInfo(list<Node>& node_list) {
	// I will create a list with elements f_i = percent of nodes with i degree
	// Initially I will make it size 37, but it can expand if necessary
	int node_tally = node_list.size();

	vector<vector<double>> graph_data(node_tally, vector<double>(node_tally));

	list<Node>::iterator itr;
	int i = 0;
	for (itr = node_list.begin(); itr != node_list.end(); itr++) {
		graph_data[i] = itr->ReturnOutWardEdges(node_list, itr);
		i++;
	}

	return graph_data;
}

// Pick which species exist in initial network
vector<int> GenerateInitialGenomes(int network_size, uniform_real_distribution<double> dist, default_random_engine gen) {
	//// Initialize vector of bools of a size equal to the number of possible species. Each bool represents whether
	//// that species has been choosen yet to exist in the initial network. Initially all bools are false.
	vector<bool> genome_options(genome_space, false);
	//// Initialize a vector of ints that are the decimal representations of the genomes of species that will exist
	//// in the initial network.
	vector<int> genome_list(network_size, 0);
	//// Initialize an int to be the random number
	int temp_rand;
	//// Do the following once for each node in the initial network
	for (int i = 0; i < network_size; i++) {
		//// ( temp_rand = (dist(gen) * genome_space) + 1) always evaluates as true because the range does not include 0
		//// (genome_options[temp_rand - 1]) only evaluates true if temp_rand has been previously choosen
		while ((temp_rand = (dist(gen) * genome_space) + 1) && (genome_options[temp_rand - 1])) {}
		//// Don't allow this genome to be picked again
		genome_options[temp_rand - 1] = true;
		//// Store the choosen decimal representation of the genome
		genome_list[i] = temp_rand - 1;
	}
	return genome_list;
}

// Allows my iterator to loop back to begining
//// loop range = 1 by default (See declaration)
void IterateCyclically(list<Node>::iterator& my_iterator, list<Node>& my_list, int loop_range) {
	for (int i = 0; i < loop_range; ++i) {
		++my_iterator;
		if (my_iterator == my_list.end()) { my_iterator = my_list.begin(); }
	}
}

// Checks if any recently updated nodes need to be deleted
// For each node that needs deletion, it calls Disentangle and then deletes it.
bool UpdateNodeList(list<Node>& node_list, bool zero_fit_death) {

	// If dead_nodes_found is true by the end of this function call, UpdateNodeList needs to be run again.
	bool dead_nodes_found = false;

	// Search for dead nodes
	//// Initialize one iterator to iterate through all nodes and one to locate the best current candidate for death
	list<Node>::iterator temp_itr = node_list.begin();
	list<Node>::iterator potential_dead_node_itr;

	for (int i = 0; i < node_list.size(); ++i) {
		//// If the node has been flagged as recently checked it may be dead
		if (temp_itr->CheckRecentlyUpdated()) {
			////// Update the node's sum
			temp_itr->UpdateSum();
			
			////// zero_fit_death comes from the commandline. If it is true, a node with 0 fitness will die if all other
			////// nodes have positive fitnesses. If it is false, only nodes with negative fitnesses can die
			if (temp_itr->ReturnSum() < 0.0 || (temp_itr->ReturnSum() == 0.0 && zero_fit_death)) {
				////// The first node that meets the conditions necessary for death is our first candidate for death
				if (!dead_nodes_found) {
					dead_nodes_found = true;
					potential_dead_node_itr = temp_itr++;
				}
				////// If a node meets the conditions necessary for death and it is not the first such node, it must be compared
				////// to the current candidate for death. If its fitness is lower than the current candidate, it becomes the current
				////// candidate. If not, we ignore it for the rest of this function call
				else {
					if (temp_itr->ReturnSum() < potential_dead_node_itr->ReturnSum()) {
						potential_dead_node_itr = temp_itr;
					}
					temp_itr++;
				}
			}
			////// If a node does not meet the conditions necessary for death it will be ignored for the remainder of this function call
			else { (temp_itr++)->ResetRecentlyUpdated(); }
		}
		//// If the node has not been flagged as recently checked it is not dead so no need to do anything with it
		else { temp_itr++; }
	}


	// If there is a dead node, delete it
	if (dead_nodes_found) {
		//// Remove pointers to other nodes
		potential_dead_node_itr->Disentangle();
		//// Remove node from node list
		node_list.erase(potential_dead_node_itr);
	}

	return dead_nodes_found;
}

// This function calculates the matrix elements of the interaction matrix
double GetCoefficient(bitset<genome_size> genome_a, bitset<genome_size> genome_b, vector<double>& X, vector<double>& Y,
	vector<bool>& connect_x, vector<bool>& connect_y) {

	// Performs the exclusive or operation between the two genome bitsets
	bitset<genome_size> genome_aXORb = genome_a;
	genome_aXORb ^= genome_b;

	// Stores the decimal representations of the two genomes and the exclusive or of them into ints
	int genome_a_dec = genome_a.to_ulong();
	int genome_b_dec = genome_b.to_ulong();
	int genome_aXORb_dec = genome_aXORb.to_ulong();

	double coefficient = 0;
	// Calculate the matrix element
	//// Define:
	////// a^b = decimal representation of exclusive or operation between genome a and genome b
	////// a = decimal representation of genome a
	////// b = decimal representation of genome b
	//// If and only if connect_x[a^b + 2(a+1)] and connect_y[a^b + 2(b+1)] are both true, there is a nonzero matrix element
	//// If there is a nonzero matrix element, it is (X[a^b + 2(a+1)] + Y[a^b + 2(a+1)])/2
	if (connect_x[genome_aXORb_dec + 2 * (genome_a_dec + 1)] && connect_y[genome_aXORb_dec + 2 * (genome_b_dec + 1)]) {
		coefficient = (X[genome_aXORb_dec + 2 * (genome_a_dec + 1)] + Y[genome_aXORb_dec + 2 * (genome_b_dec + 1)]) / 2;
	}

	return coefficient;
}

// I didn't comment out this function in depth, but it essentially calculates the correlation between the coefficients of
// interaction for pairs of nodes as a function of the Hamming distance between the pairs
void CheckCorrelations(vector<double>& x, vector<double>& y, vector<bool>& x_connect, vector<bool>& y_connect) {
	vector<double> same(2 * genome_size + 1, 0);
	vector<double> coeff_prod(2 * genome_size + 1, 0);
	vector<double> mean_coeff_1(2 * genome_size + 1, 0);
	vector<double> mean_value_1(2 * genome_size + 1, 0);
	vector<double> coeff_variance(2 * genome_size + 1, 0);
	vector<double> connect_variance(2 * genome_size + 1, 0);
	double mean_coeff = 0;
	double mean_value = 0;
	double total_amount = 0;
	vector<double> total(2 * genome_size + 1, 0);
	for (int i = 0; i < genome_space; ++i) {
		for (int j = 0; j < genome_space; ++j) {
			if (i != j) {
				bitset<genome_size> one_genome = i;
				bitset<genome_size> other_genome = j;
				bitset<genome_size> mixed_genome = one_genome;
				mixed_genome ^= other_genome;
				int number = mixed_genome.to_ulong();
				bool connect;
				double coeff_1 = (x[number + 2 * (i + 1)] + y[number + 2 * (j + 1)]) / 2;

				total_amount += 1;
				if (x_connect[number + 2 * (i + 1)] && y_connect[number + 2 * (j + 1)]) {
					connect = true;
				}
				else {
					connect = false;
				}


				for (int l = 0; l < genome_space; ++l) {
					for (int m = 0; m < genome_space; ++m) {
						if ((l != m) && (genome_space * i + j >= genome_space * l + m)) {
							bitset<genome_size> one_genome_2 = l;
							bitset<genome_size> other_genome_2 = m;
							int number_b = other_genome_2.to_ulong();
							bitset<genome_size> mixed_genome = one_genome_2;
							mixed_genome ^= other_genome_2;
							int number = mixed_genome.to_ulong();
							bool connect_2;
							double coeff_2 = (x[number + 2 * (l + 1)] + y[number + 2 * (m + 1)]) / 2;


							if (x_connect[number + 2 * (l + 1)] && y_connect[number + 2 * (m + 1)]) {
								connect_2 = true;
							}
							else { connect_2 = false; }
							int h_dist = (one_genome^one_genome_2).count() + (other_genome^other_genome_2).count();
							total[h_dist] += 1;
							mean_coeff_1[h_dist] += coeff_1 + coeff_2;

							if (connect) { mean_value_1[h_dist] += 1; }
							else { mean_value_1[h_dist] -= 1; }

							if (connect_2) { mean_value_1[h_dist] += 1; }
							else { mean_value_1[h_dist] -= 1; }

							connect_variance[h_dist] += (2 * connect - 1) * (2 * connect - 1);

							coeff_prod[h_dist] += coeff_1 * coeff_2;
							if (connect == connect_2) { same[h_dist] += 1; }
							else { same[h_dist] -= 1; }

							coeff_variance[h_dist] += coeff_1 * coeff_1;
						}
					}
				}
			}
		}
	}


	for (int i = 0; i < 2 * genome_size + 1; i++) {
		mean_value_1[i] /= (2 * total[i]);
		same[i] /= total[i];
		connect_variance[i] = connect_variance[i] / total[i] - mean_value_1[i] * mean_value_1[i];

		mean_coeff_1[i] /= (2 * total[i]);
		coeff_prod[i] /= total[i];
		coeff_variance[i] = coeff_variance[i] / total[i] - mean_coeff_1[i] * mean_coeff_1[i];
	}

	ofstream coeff_correlations_file("coefficient_correlations.txt");
	//coeff_correlations_file << "mean value of correlations: " << mean_coeff << "\n";

	ofstream connect_correlations_file("connectedness_correlations.txt");
	//connect_correlations_file << "mean value of connectedness pairs: " << mean_value << "\n";

	for (int i = 0; i < 2 * genome_size + 1; i++) {
		connect_correlations_file << i << "  " << total[i] << " " << (same[i] - mean_value_1[i] * mean_value_1[i]) / connect_variance[i] << "\n";
	}

	for (int i = 0; i < 2 * genome_size + 1; i++) {
		coeff_correlations_file << i << "  " << total[i] << " " << (coeff_prod[i] - mean_coeff_1[i] * mean_coeff_1[i]) / coeff_variance[i] << "\n";
	}

	connect_correlations_file.close();
	coeff_correlations_file.close();
}

#include "Functions.h"
#include <set>

using namespace std;

void PrintOutMessage(int message_id) {
	if (message_id == 1) { cerr << "generated nodes of initial community" << endl; }
	else if (message_id == 2) { cerr << "assigned genomes to initial community" << endl; }
	else if (message_id == 3) { cerr << "added edges to initial community" << endl; }
	else if (message_id == 4) { cerr << "updated initial community" << endl; }
}

// Pick which species exist in initial network
vector<int> GenerateInitialGenomes(int network_size, uniform_real_distribution<double> dist, default_random_engine gen, long genome_space) {
	set<int> s;
	vector<int> genome_list;
	while( s.size() < network_size ) {
		int r = static_cast<int>(genome_space * dist(gen));
    if( s.find(r) == s.end() ) {
      genome_list.push_back(r);
			s.insert(r);
		}
	}
	return std::move(genome_list);
}

// Checks if any recently updated nodes need to be deleted
// For each node that needs deletion, it calls Disentangle and then deletes it.
bool UpdateNodeList(list<Node>& node_list, bool zero_fit_death, long step) {
	// Search for dead nodes
	//// Initialize one iterator to iterate through all nodes and one to locate the best current candidate for death
	auto smallest = node_list.begin();
  for(auto it = node_list.begin(); it != node_list.end(); ++it) {
    if( it->CheckRecentlyUpdated() ) {
			it->UpdateSum();
      double f = it->ReturnSum();
			if( f > 0 ) { it->ResetRecentlyUpdated(); } // the species doesn't have a negative fitness
			if(f < smallest->ReturnSum()) { smallest = it; };
		}
	}

	if(smallest->ReturnSum() < 0.0 || (smallest->ReturnSum() == 0.0 && zero_fit_death)) {
		//// If past burn-in period, record lifespan of dead node
		smallest->UpdateLifespanDistribution(step);

		//// Remove pointers to other nodes
		smallest->Disentangle();
		//// Remove node from node list
		node_list.erase(smallest);
		return true;
	}

	return false;
}

// This function calculates the matrix elements of the interaction matrix
double GetCoefficient(const genome_t a, const genome_t b, const vector<double>& X, const vector<double>& Y,
	const vector<bool>& connect_x, const vector<bool>& connect_y) {

	// Performs the exclusive or operation between the two genome bitsets
	genome_t aXORb = a ^ b;

	// Stores the decimal representations of the two genomes and the exclusive or of them into ints
	unsigned long a_dec = a.to_ulong();
	unsigned long b_dec = b.to_ulong();
	unsigned long ab_dec = aXORb.to_ulong();

	double coefficient = 0.0;
	// Calculate the matrix element
	//// Define:
	////// a^b = decimal representation of exclusive or operation between genome a and genome b
	////// a = decimal representation of genome a
	////// b = decimal representation of genome b
	//// If and only if connect_x[a^b + 2(a+1)] and connect_y[a^b + 2(b+1)] are both true, there is a nonzero matrix element
	//// If there is a nonzero matrix element, it is (X[a^b + 2(a+1)] + Y[a^b + 2(a+1)])/2
	if (connect_x[ab_dec + 2 * (a_dec + 1)] && connect_y[ab_dec + 2 * (b_dec + 1)]) {
		coefficient = (X[ab_dec + 2 * (a_dec + 1)] + Y[ab_dec + 2 * (b_dec + 1)]) / 2;
	}

	return coefficient;
}

// I didn't comment out this function in depth, but it essentially calculates the correlation between the coefficients of
// interaction for pairs of nodes as a function of the Hamming distance between the pairs
void CheckCorrelations(vector<double>& x, vector<double>& y, vector<bool>& x_connect, vector<bool>& y_connect, long genome_size) {
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
  long genome_space = 1 << genome_size;
	for (int i = 0; i < genome_space; ++i) {
		for (int j = 0; j < genome_space; ++j) {
			if (i != j) {
				genome_t one_genome = i;
				genome_t other_genome = j;
				genome_t mixed_genome = one_genome;
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
							genome_t one_genome_2 = l;
							genome_t other_genome_2 = m;
							int number_b = other_genome_2.to_ulong();
							genome_t mixed_genome = one_genome_2;
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

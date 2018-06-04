#pragma once

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "Node.h"

using namespace std;

void PrintOutMessage(int message_id);

// Pick which species exist in initial network
vector<int> GenerateInitialGenomes(int network_size, uniform_real_distribution<double> dist, default_random_engine gen, long genome_space);

// Checks if any recently updated nodes need to be deleted
// For each node that needs deletion, it calls Disentangle and then deletes it.
bool UpdateNodeList(list<Node>& node_list, bool zero_fit_death, long step);

// This function calculates the matrix elements of the interaction matrix
double GetCoefficient(const genome_t a, const genome_t b, const vector<double>& X, const vector<double>& Y,
	const vector<bool>& connect_x, const vector<bool>& connect_y);

// I didn't comment out this function in depth, but it essentially calculates the correlation between the coefficients of
// interaction for pairs of nodes as a function of the Hamming distance between the pairs
void CheckCorrelations(vector<double>& x, vector<double>& y, vector<bool>& x_connect, vector<bool>& y_connect, long genome_size);

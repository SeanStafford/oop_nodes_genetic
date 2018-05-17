#pragma once

#include "CompileTimeConstants.h"

using namespace std;


double GetCoefficient(bitset<genome_size>, bitset<genome_size>, vector<double>&, vector<double>&, vector<bool>&, vector<bool>&);

void IterateCyclically(list<Node>::iterator&, list<Node>&, int loop_range = 1);

bool UpdateNodeList(list<Node>&, bool);

void GetModelParameters(double(&double_array)[num_input_params], int, char**);

void PrintOutMessage(int);

string AssignFileName(double(&double_array)[num_input_params]);

vector<vector<float>> MeasureDegreeDistribution(list<Node>&);

vector<vector<double>> RetrieveGraphInfo(list<Node>&);

void CheckCorrelations(vector<double>&, vector<double>&, vector<bool>&, vector<bool>&);

vector<int> GenerateInitialGenomes(int, uniform_real_distribution<double>, default_random_engine);
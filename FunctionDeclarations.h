#pragma once

using namespace std;


double GetCoefficient(genome_t, genome_t, vector<double>&, vector<double>&, vector<bool>&, vector<bool>&);

void IterateCyclically(list<Node>::iterator&, list<Node>&, int loop_range = 1);

bool UpdateNodeList(list<Node>&, bool);

void PrintOutMessage(int);

vector<vector<float>> MeasureDegreeDistribution(list<Node>&);

vector<vector<double>> RetrieveGraphInfo(list<Node>&);

void CheckCorrelations(vector<double>&, vector<double>&, vector<bool>&, vector<bool>&);

vector<int> GenerateInitialGenomes(int, uniform_real_distribution<double>, default_random_engine);
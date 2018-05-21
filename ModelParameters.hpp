//
// Created by Yohsuke Murase on 2018/05/21.
//

#ifndef OOP_NODES_GENETIC_MODELPARAMETERS_HPP
#define OOP_NODES_GENETIC_MODELPARAMETERS_HPP

#include <iostream>

class ModelParameters {
public:
  ModelParameters(bool check_cor, bool deg_fixed, bool zero_f_ext, bool edge_direct, long init_pop_size, long deg,
                  double density, long total_steps, long output_steps, long s) :
      checking_correlation(check_cor),
      degree_fixed(deg_fixed),
      zero_fitness_extinction(zero_f_ext),
      edge_direction(edge_direct),
      initial_population_size(init_pop_size),
      degree(deg),
      link_density(density),
      total_steps(total_steps),
      output_steps(output_steps),
      seed(s) {};
  const bool checking_correlation;
  const bool degree_fixed;
  const bool zero_fitness_extinction;  // true: make zero-fitness species go extinct. false: let it survive.
  const bool edge_direction; // [TODO] not implemented?
  const long initial_population_size;
  const long degree;      // valid only if fixed_degree is true
  const double link_density; // valid only if fixed_degree is false
  const long total_steps;
  const long output_steps;
  const long seed;

  void Print() {
    std::cout
        << "checking_correlation : " << checking_correlation << std::endl
        << "degree_fixed : " << degree_fixed << std::endl
        << "zero_fitness_extinction : " << zero_fitness_extinction << std::endl
        << "edge_direction : " << edge_direction << std::endl
        << "initial_population_size : " << initial_population_size << std::endl
        << "degree : " << degree << std::endl
        << "link_density : " << link_density << std::endl
        << "total_steps : " << total_steps << std::endl
        << "output_steps : " << output_steps << std::endl;
  }

  static ModelParameters Load(int argc, char** argv) {
    if( argc != 10 ) {
      std::cerr << "[Error] invalid number of arguments" << std::endl;
      std::cerr << "  Usage : ./a.out "
          << "<1:checking_correlation:bool> "
          << "<2:degree_fixed_or_density:bool> "
          << "<3:zero_fitness_extinction:bool> "
          << "<4:edge_direction:bool> "
          << "<5:initial_population_size:int> "
          << "<6:fixed_degree:int or degree_density:double> "
          << "<7:total_steps:int> "
          << "<8:output_step:int> "
          << "<9:seed:int> "
          << std::endl;
    }
    bool check_corr = atoi(argv[1]) != 0;
    bool degree_fixed = atoi(argv[2]) != 0;
    bool zero_fitness_ext = atoi(argv[3]) != 0;
    bool edge_direct = atoi(argv[4]) != 0;
    long init_pop_size = atol(argv[5]);
    long k = 0;
    double c = 0.0;
    if( degree_fixed ) {
      k = atol(argv[6]);
    }
    else {
      c = atof(argv[6]);
    }
    long total_steps = atol(argv[7]);
    long output_steps = atol(argv[8]);
    long seed = atol(argv[9]);

    ModelParameters param(check_corr, degree_fixed, zero_fitness_ext, edge_direct, init_pop_size, k, c, total_steps, output_steps, seed);
    return std::move(param);
  }
};

#endif //OOP_NODES_GENETIC_MODELPARAMETERS_HPP

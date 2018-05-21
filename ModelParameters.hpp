//
// Created by Yohsuke Murase on 2018/05/21.
//

#ifndef OOP_NODES_GENETIC_MODELPARAMETERS_HPP
#define OOP_NODES_GENETIC_MODELPARAMETERS_HPP

#include <iostream>

class ModelParameters {
public:
  ModelParameters(bool check_cor, bool zero_f_ext, long init_pop_size,
                  double density, long total_steps, long output_steps, long s) :
      checking_correlation(check_cor),
      zero_fitness_extinction(zero_f_ext),
      initial_population_size(init_pop_size),
      link_density(density),
      total_steps(total_steps),
      output_steps(output_steps),
      seed(s) {};
  const bool checking_correlation;
  const bool zero_fitness_extinction;  // true: make zero-fitness species go extinct. false: let it survive.
  const long initial_population_size;
  const double link_density; // valid only if fixed_degree is false
  const long total_steps;
  const long output_steps;
  const long seed;

  void Print() {
    std::cout
        << "checking_correlation : " << checking_correlation << std::endl
        << "zero_fitness_extinction : " << zero_fitness_extinction << std::endl
        << "initial_population_size : " << initial_population_size << std::endl
        << "link_density : " << link_density << std::endl
        << "total_steps : " << total_steps << std::endl
        << "output_steps : " << output_steps << std::endl;
  }

  static ModelParameters Load(int argc, char** argv) {
    if( argc != 8 ) {
      std::cerr << "[Error] invalid number of arguments" << std::endl;
      std::cerr << "  Usage : ./a.out "
          << "<1:checking_correlation:bool> "
          << "<2:zero_fitness_extinction:bool> "
          << "<3:initial_population_size:int> "
          << "<4:degree_density:double> "
          << "<5:total_steps:int> "
          << "<6:output_step:int> "
          << "<7:seed:int> "
          << std::endl;
    }
    bool check_corr = atoi(argv[1]) != 0;
    bool zero_fitness_ext = atoi(argv[2]) != 0;
    long init_pop_size = atol(argv[3]);
    double c = atof(argv[4]);
    long total_steps = atol(argv[5]);
    long output_steps = atol(argv[6]);
    long seed = atol(argv[7]);

    ModelParameters param(check_corr, zero_fitness_ext, init_pop_size, c, total_steps, output_steps, seed);
    return std::move(param);
  }
};

#endif //OOP_NODES_GENETIC_MODELPARAMETERS_HPP

//
// Created by Yohsuke Murase on 2018/05/21.
//

#ifndef OOP_NODES_GENETIC_MODELPARAMETERS_HPP
#define OOP_NODES_GENETIC_MODELPARAMETERS_HPP

#include <iostream>

class ModelParameters {
public:
  ModelParameters(bool check_cor, bool zero_f_ext, long genome_l, long init_pop_size,
                  double density, long total_steps, long output_steps, double mean, long burn_in, long s) :
      checking_correlation(check_cor),
      zero_fitness_extinction(zero_f_ext),
      genome_length(genome_l),
      initial_population_size(init_pop_size),
      link_density(density),
      total_steps(total_steps),
      output_steps(output_steps),
	  mean(mean),
	  burn_in(burn_in),
      seed(s) {};
  const bool checking_correlation;
  const bool zero_fitness_extinction;  // true: make zero-fitness species go extinct. false: let it survive.
  const long genome_length;
  const long initial_population_size;
  const double link_density; // valid only if fixed_degree is false
  const long total_steps;
  const long output_steps;
  const double mean; // mean of distribution from which entries in X and Y are drawn
  const long burn_in; // after how many time steps do we start keeping track of lifespans
  const long seed;
  long GenomeSpace() const {
    return 1 << genome_length;
  }

  void Print() {
	  std::cerr
		  << "checking_correlation : " << checking_correlation << std::endl
		  << "zero_fitness_extinction : " << zero_fitness_extinction << std::endl
		  << "genome_length : " << genome_length << std::endl
		  << "initial_population_size : " << initial_population_size << std::endl
		  << "link_density : " << link_density << std::endl
		  << "total_steps : " << total_steps << std::endl
		  << "output_steps : " << output_steps << std::endl
		  << "coefficient_distribution_mean : " << mean << std::endl
		  << "lifespan_burn_in : " << burn_in << std::endl;

  }

  static ModelParameters Load(int argc, char** argv) {
    if( argc != 11 ) {
      std::cerr << "[Error] invalid number of arguments" << std::endl;
      std::cerr << "  Usage : ./a.out "
          << "<1:checking_correlation:bool> "
          << "<2:zero_fitness_extinction:bool> "
          << "<3:genome_length:int> "
          << "<4:initial_population_size:int> "
          << "<5:degree_density:double> "
          << "<6:total_steps:int> "
          << "<7:output_step:int> "
		  << "<8:coefficient_distribution_mean:double> "
		  << "<9:lifespan_burn_in:int> "
          << "<10:seed:int> "
          << std::endl;
      throw "invalid number of arguments";
    }
    bool check_corr = atoi(argv[1]) != 0;
    bool zero_fitness_ext = atoi(argv[2]) != 0;
    long genome_length = atoi(argv[3]);
    long init_pop_size = atol(argv[4]);
    double c = atof(argv[5]);
    long total_steps = atol(argv[6]);
    long output_steps = atol(argv[7]);
	double mean = atof(argv[8]);
	long burn_in = atol(argv[9]);
    long seed = atol(argv[10]);

    if( init_pop_size > (1<<genome_length)) {
      throw "init_pop_size is too large";
    }
    if( total_steps < output_steps ) {
      throw "output_steps must be smaller than total_steps";
    }
    if( genome_length > 64) {
      throw "genome length cannot be longer than 64";
    }
    if( c < 0.0 || c > 1.0) {
      throw "c must be in [0,1]";
    }

    ModelParameters param(check_corr, zero_fitness_ext, genome_length, init_pop_size, c, total_steps, output_steps, mean, burn_in, seed);
    return std::move(param);
  }
};

#endif //OOP_NODES_GENETIC_MODELPARAMETERS_HPP

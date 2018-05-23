//
// Created by Yohsuke Murase on 2018/05/23.
//

#ifndef OOP_NODES_GENETIC_RESULTS_HPP
#define OOP_NODES_GENETIC_RESULTS_HPP

#include <iostream>
#include <fstream>

struct Results {
  long initial_population_size;
  double average_population_size;
  long final_population_size;
  long total_steps; // steps until the whole extinction
  void PrintJson(std::ostream& out) {
    out << "{\n"
        << "\"initial_population_size\": " << initial_population_size << ",\n"
        << "\"average_population_size\": " << average_population_size << ",\n"
        << "\"final_population_size\": " << final_population_size << ",\n"
        << "\"total_steps\": " << total_steps << "\n"
        << "}";
    out.flush();
  }
};

#endif //OOP_NODES_GENETIC_RESULTS_HPP

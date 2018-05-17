#pragma once

// It may be desirable to save this information in a JSON document and have the OACIS simulator compile and then run my program everytime.
// That way I could run different genome sizes from the same simulator.
// However I don't know how to do that and having to compile everytime may cause the simulator to take longer than necessary
const int genome_size = 8;
const int num_input_params = 9;
const int genome_space = pow(2, genome_size);
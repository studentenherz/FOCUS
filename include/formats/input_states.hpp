#if !defined(FOCUS_INCLUDE_FORMATS_INPUT_STATES_HPP)
#define FOCUS_INCLUDE_FORMATS_INPUT_STATES_HPP

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "regex_tokenizer.hpp"
#include "types/array.hpp"
#include "types/vector.hpp"


/**
 * Read initial states of particles from a file in the format
 * r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s]
 *
 * Comments start with #
 *  
 * @return Array of initial states
 */
Array<State> load_states(std::string filename){
	std::ifstream fi(filename);
	if(!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return Array<State>(0);
	}

	Tokenizer<std::ifstream> tk("[+-]?\\d*[\\.]?\\d+(?:[Ee][+-]?\\d+)?"); // captures any number;
	std::string token;

	std::vector<State> v_states;

	double r, q, z, vr, vq, vz;
	std::string line;
	while	(std::getline(fi, line)){
		if (regex_match(line, "#.*") || line.empty()) continue;
		std::istringstream iss(line);
		iss >> r >> q >> z >> vr >> vq >> vz;
		v_states.push_back({r, q, z, vr, vq, vz});
	}

	Array<State> states(v_states.size());
	for (size_t i = 0; i < states.size(); i++)
		states[i] = v_states[i];

	return states;
}

#endif // FOCUS_INCLUDE_FORMATS_INPUT_STATES_HPP

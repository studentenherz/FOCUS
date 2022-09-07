#if !defined(FOCUS_INCLUDE_FORMATS_PARTICLE_STATES_HPP)
#define FOCUS_INCLUDE_FORMATS_PARTICLE_STATES_HPP

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "regex_tokenizer.hpp"
#include "types/array.hpp"
#include "types/vector.hpp"
#include "types/particle.hpp"

/**
 * Read initial states of particles from a file in the format
 * r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s]
 *
 * Comments start with #
 *  
 * @param filename input file
 * @return Array of initial states
 */
Array<State> load_states(std::string filename){
	std::ifstream fi(filename);
	if(!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return Array<State>(0);
	}

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

/**
 * Dump states of particles into a file in the format
 * r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s]
 *
 * Comments start with #
 *  
 * @param filename output file
 * @param states of the particles
 * @param part particle species
 * @return true on success
 */
bool dump_states(std::string filename, Array<State>& states, Particle part){
	std::ofstream fo(filename);
	if(!fo.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	#ifdef GIT_VERSION
	fo << "# Output file created with FOCUS git version of the code " << GIT_VERSION << '\n';
	#else
	fo << "# Output file created with FOCUS git version of the code unknown\n";
	#endif
	fo << "# The format of the output is as follows:\n";
	fo << "# r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s]\n\n";
	fo << "# Particle species:\n";
	fo << "# Z = " << part.q << '\n';
	fo << "# m = " << part.m << " Da\n\n";
	fo << "#States:\n";

	for (size_t i = 0; i < states.size(); i++)
		fo << states[i] << '\n';

	fo.close();
	return true;
}

#endif // FOCUS_INCLUDE_FORMATS_PARTICLE_STATES_HPP

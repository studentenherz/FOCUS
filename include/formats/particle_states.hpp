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
		if (regex_match(line, "#.*") || line.empty() || line.find("nan") != line.npos) continue;
		std::istringstream iss(line);
		iss >> r >> q >> z >> vr >> vq >> vz;
		State stt{r, q, z, vr, vq, vz};
		v_states.push_back(stt);
	}

	Array<State> states(v_states.size());
	for (size_t i = 0; i < states.size(); i++)
		states[i] = v_states[i];

	return states;
}

/**
 * Read initial states of particles from a file in the format
 * r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s]
 *
 * Comments start with #
 *  
 * @param filename input file
 * @return Array of initial states
 */
void load_states(std::string filename, Array<State>& states){
	std::ifstream fi(filename);
	if(!fi.is_open()){
		std::cerr << "\x1b[31mError\x1b[0m opening file \x1b[1m" << filename << "\x1b[0m!\n";
		return;
	}

	std::vector<State> v_states;

	double r, q, z, vr, vq, vz;
	std::string line;
	while	(std::getline(fi, line)){
		if (regex_match(line, "#.*") || line.empty() || line.find("nan") != line.npos) continue;
		std::istringstream iss(line);
		iss >> r >> q >> z >> vr >> vq >> vz;
		State stt{r, q, z, vr, vq, vz};
		v_states.push_back(stt);
	}

	states.resize(v_states.size());
	for (size_t i = 0; i < states.size(); i++)
		states[i] = v_states[i];
}

template<typename T>
void print_ith_line(size_t i, std::ostream& os, T arr){
	os << arr[i] << ' ';
}

template<typename T, typename... Args>
void print_ith_line(size_t i, std::ostream& os, T arr, Args... args){
	os << arr[i] << ' ';
	print_ith_line(i, os, args...);
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
template<typename... Args>
bool dump_states(std::string filename, Array<State>& states, Particle part, std::string desc = "", Args&... args){
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
	fo << "# r [m] | theta [rad] | z [m] | vr [m/s] | vtheta [m/s] | vz [m/s] | " << desc << "\n\n";
	fo << "# Particle species:\n";
	fo << "# Z = " << part.q << '\n';
	fo << "# m = " << part.m << " Da\n\n";
	fo << "#States:\n";

	for (size_t i = 0; i < states.size(); i++){
		print_ith_line(i, fo, states, args...);
		fo << '\n';
	}

	fo.close();
	return true;
}

#endif // FOCUS_INCLUDE_FORMATS_PARTICLE_STATES_HPP

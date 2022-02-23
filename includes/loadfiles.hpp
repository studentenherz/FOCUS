#if !defined(FOCUS_LOADFILES_HPP)
#define FOCUS_LOADFILES_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "types.hpp"

/**
 * Load a matrix from a file expected to have a row on each line.
 * @param filename The name of the file.
 * @return true if the file was read successfully, false otherwise.
 */
bool load_matrix_from_file(std::string filename, Matrix2D& M){
	std::ifstream fi(filename);
	if (!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	std::string line;
	// for each line
	while(getline(fi, line)){
		std::istringstream ss(line);
		M.push_back({});
		double value;
		// for each value in line
		while(ss >> value) M[M.size() - 1].push_back(value);
	}
	
	if (M.empty()) std::cerr << "File " << filename << " is empty\n";
	return !M.empty();
}

#endif // FOCUS_LOADFILES_HPP

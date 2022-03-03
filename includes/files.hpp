#if !defined(FOCUS_FILES_HPP)
#define FOCUS_FILES_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "types.hpp"

/**
 * Load a matrix from a file expected to have a row on each line.
 * @param filename The name of the file.
 * @param M the matrix to store the data.
 * @return true if the file was read successfully, false otherwise.
 */
bool load(std::string filename, Matrix2D& M){
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

/**
 * Load a matrix of dimensions (Nx, Ny) from file.
 * @param filename The name of the file.
 * @param M the matrix to store the data.
 * @param Nx dimension in x.
 * @param Ny dimension in y.
 * @return true if the file was read successfully, false otherwise.
 */
bool load(std::string filename, Matrix2D& M, int Nx, int Ny){
	std::ifstream fi(filename);
	if (!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	double v;
	for (int i = 0; i < Nx; i++){
		M.push_back({});
		for (int j = 0; j < Ny; j++){
			fi >> v;
			M[i].push_back(v);
		}
	}
	
	if (M.empty()) std::cerr << "File " << filename << " is empty\n";
	return !M.empty();
}

/**
 * Load a matrix from a file expected to have a row on each line.
 * @param filename The name of the file.
 * @return true if the file was read successfully, false otherwise.
 */
bool dump(std::string filename, Matrix2D& M){
	std::ofstream fo(filename);
	if (!fo.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	for(auto line : M){
		for (auto v : line){
			fo << v << ' ';
		}
		fo << '\n';
	}
	return true;
}

#endif // FOCUS_FILES_HPP

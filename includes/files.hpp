#if !defined(FOCUS_FILES_HPP)
#define FOCUS_FILES_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "types/matrix2D.hpp"

/**
 * Load a matrix of dimensions (Nx, Ny) from file. The file is 
 * expected to have the dimensions of the matrix as the first 
 * two integers.
 * @param filename the name of the file.
 * @param M the matrix to store the data.
 * @param Nx dimension in x.
 * @param Ny dimension in y.
 * @return true if the file was read successfully, false otherwise.
 */
template <typename T>
bool load(std::string filename, Matrix2D<T>& M){
	std::ifstream fi(filename);
	if (!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	size_t Nx, Ny;
	fi >> Nx >> Ny;

	M.reshape(Nx, Ny);

	double v;
	for (size_t j = 0; j < Ny; j++){
		for (size_t i = 0; i < Nx; i++){
			if(!(fi >> v)) return false;
			M(i, j) = v;
		}
	}
	
	return true;
}

/**
 * Dump the content of a matrix to a file.
 * @param filename the name of the file.
 * @param M the matrix.
 * @param matrix_shape if true the shape of the matrix will be written too.
 * @return true if the file was written successfully, false otherwise.
 */
template<typename T>
bool dump(std::string filename, const Matrix2D<T>& M, bool matrix_shape = true){
	std::ofstream fo(filename);
	if (!fo.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return false;
	}

	shape_t shape = M.shape();
	if (matrix_shape)
		fo << shape.first << ' ' << shape.second << '\n';

	for (size_t j = 0; j < shape.second; j++){
		for(size_t i = 0; i < shape.first; i++){
			fo << M(i, j) << ' ';
		}
		fo << '\n';
	}
	return true;
}

#endif // FOCUS_FILES_HPP

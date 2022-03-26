#include <iostream>
#include <iomanip>

#include "Chebyshev.hpp"
#include "files.hpp"

int main(int argc, char const *argv[]){
	if (argc < 4){
		std::cout << "Usage:\nchebyshev_expansion <file_in> <order> <error>\n";
		return 0;
	}

	// Raw EFIT matrix
	Matrix2D<double> M;
	if (!load(argv[1], M))
		return 1;

	// Chebyshev coefficients vector
	size_t n = atoi(argv[2]);
	Matrix2D<double> a(n + 1, n + 1);

	// Just to set some values,
	// these are unimportant right now
	const double x_min = 0;
	const double x_max = 1;
	const double y_min = 0;
	const double y_max = 1;

	// Get expansion coefficients
	Chebyshev_T_expansion(n, a, M, x_min, x_max, y_min, y_max);	

	// Calculate same Matrix from expansion
	size_t Nx, Ny;
	Nx = M.shape().first;
	Ny = M.shape().second;
	Matrix2D<double> new_M(Nx, Ny);

	for (size_t i = 0; i < Nx; i++){
		double x = x_min + i * (x_max - x_min) / Nx;
		for (size_t j = 0; j < Ny; j++){
			double y = y_min + j * (y_max - y_min) / Ny;
			new_M(i, j) = evaluate_Chebyshev_T_expansion(n, a, x, y, x_min, x_max, y_min, y_max);
		}
	}

	// Compare
	double epsilon = atof(argv[3]);
	double max_diff = 0;
	for (size_t i = 0; i < Nx; i++)
		for (size_t j = 0; j < Ny; j++){
			max_diff = std::max(max_diff, std::abs(M(i, j) - new_M(i, j)));
			if (max_diff > epsilon){
				std::cerr << "At position (" << i << ", " << j << ") got " << new_M(i, j) << " but expected " << M(i, j) << '\n';
				return 1;
			}
		}

	std::cout << "Max difference: " << max_diff << '\n';

	return 0;
}

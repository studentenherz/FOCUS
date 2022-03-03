#include <iostream>
#include <iomanip>

#include "Chebyshev.hpp"
#include "files.hpp"

int main(int argc, char const *argv[]){
	if (argc < 3) return 0;

	// Raw EFIT matrix
	Matrix2D M;
	if (!load(argv[1], M))
		return 1;

	// Chebyshev coefficients vector
	int n = atoi(argv[2]);
	Matrix2D a;
	for(int i=0; i<=n; i++)
		a.push_back(std::vector<double>(n + 1));

	// Just to set some values,
	// these are unimportant right now
	const double x_min = 0;
	const double x_max = 1;
	const double y_min = 0;
	const double y_max = 1;

	// Get expansion coefficients
	Chebyshev_T_expansion(n, a, M, x_min, x_max, y_min, y_max);	

	// Calculate same Matrix from expansion
	Matrix2D new_M;
	int Nx = M.size();
	int Ny = M[0].size();

	for(int i=0; i<Nx; i++)
		new_M.push_back(std::vector<double>(Ny));

	for (int i = 0; i < Nx; i++){
		double x = x_min + i * (x_max - x_min) / Nx;
		for (int j = 0; j < Ny; j++){
			double y = y_min + j * (y_max - y_min) / Ny;
			new_M[i][j] = evaluate_Chebyshev_T_expansion(a, x, y, x_min, x_max, y_min, y_max);
		}
	}

	dump("from_Chebyshev.out", new_M);

	return 0;
}

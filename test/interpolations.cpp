#include <iostream>

#include "files.hpp"
#include "interpolations.hpp"

int main(int argc, char const *argv[]){
	if (argc < 2) return 0;
	Matrix2D psi;
	if (load(argv[1], psi))
		std::cout << psi.size() << " x " << psi[0].size() << " matrix read\n";
	
	// gives some value
	double v = six_point_formula(0.5, 0.6, psi, 0, 1, 0, 1);
	std::cout << "psi(0.5, 0.6) = " << v << '\n';
	
	// should give nan and display error message
	v = six_point_formula(1.5, 0.6, psi, 0, 1, 0, 1);
	std::cout << "psi(1.5, 0.6) = " << v << '\n';

	return 0;
}
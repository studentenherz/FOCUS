#include <iostream>

#include "formats/matrix.hpp"
#include "interpolations.hpp"

int main(int argc, char const *argv[]){
	if (argc < 2) return 0;
	Matrix2D<double> psi;
	if (!load(argv[1], psi))
		return 1;	
	
	ScalarField f(psi, 0, 1, 0, 1);

	// gives some value
	double v = six_point_formula(0.23, 0.8, f);
	if (std::abs(v - 0.242005) > 1e-6){
		std::cerr << "Interpolation vale gave " << v << " but 0.242005 was expected (different by at least 1e-6)\n";
		return 1;
	}
	
	// should give nan
	v = six_point_formula(1.5, 0.6, f);
	if (!std::isnan(v)){
		std::cerr << "Interpolation vale gave " << v << " but nan was expected\n";
		return 1;
	}

	return 0;
}
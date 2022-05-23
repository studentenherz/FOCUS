#include <iostream>
#include <cmath>

#include "chebyshev.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "files.hpp"

double f(double x, double y){
	return 3 * sin(x / 10) * cos(2 * y + 1);
}

double df_dx(double x, double y){
	return 0.3 * cos(x / 10) * cos(2 * y + 1);
}

double df_dy(double x, double y){
	return -6 * sin(x / 10) * sin(2 * y + 1);
}

int main(int argc, char const *argv[]){
	if (argc < 2){
		std::cout << "Usage:\nchebyshev... <input_file>\n";
		return 0;
	}

	double x_min = 0.9682511999999;
	double x_max = 2.4097488000000;
	double y_min = -1.442599200000;
	double y_max = 1.442599200000;

	double rlow = 0.945960;
	double rhigh = 2.4458;
	double zlow = -1.47204;
	double zhigh = 1.47204;

	Matrix2D<double> M;
	load(argv[1], M);
	ScalarField f(M, rlow, rhigh, zlow, zhigh);

	size_t n = 26;

	ChebyshevExpansion ch(n, f, x_min, x_max, y_min, y_max);

	// Calculate the coefficients of the expansion
	size_t Nn = 600;
	Matrix2D<double> Br(Nn, Nn);
	Matrix2D<double> Bz(Nn, Nn);
	Matrix2D<double> Psi(Nn, Nn);


	// double epsilon = 0.1;

	const double B_0 = 1.77048;
	const double a_m = 1.00000;

	for(size_t i = 0; i<Nn; i++){
		double x = x_min + (x_max - x_min) * i / (Nn - 1);
		for(size_t j = 0; j<Nn; j++){
			double y = y_min + (y_max - y_min) * j / (Nn - 1);
			Br(i, j) =  - ch.dy(x, y) / x / B_0 /a_m / a_m;
			if(std::isnan(Br(i, j)))
				Br(i, j) = 0;
			Bz(i, j) =  ch.dx(x, y) / x / B_0 /a_m / a_m;
			if(std::isnan(Bz(i, j)))
				Bz(i, j) = 0;
			Psi(i, j) = ch(x, y);
		}
	}

	dump("Br.dat", Br, false);
	dump("Bz.dat", Bz, false);
	dump("Psi.dat", Psi, false);

	return 0;
}

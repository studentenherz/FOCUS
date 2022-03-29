#include <iostream>
#include <cmath>

#include "Chebyshev.hpp"
#include "types/matrix2D.hpp"
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
	if (argc < 3){
		std::cout << "Usage:\nchebyshev... <input_file> <output_file>\n";
		return 0;
	}

	double x_min = 0.9682511999999;
	double x_max = 2.4097488000000;
	double y_min = -1.442599200000;
	double y_max = 1.442599200000;

	Matrix2D<double> M;
	load(argv[1], M);


	size_t n = 25;
	Matrix2D<double> a(n, n);

	// Calculate the coefficients of the expansion
	Chebyshev_T_expansion(n, a, M, x_min, x_max, y_min, y_max);
	size_t Nn = 400;
	Matrix2D<double> Br(Nn, Nn);


	// double epsilon = 0.1;



	for(size_t i = 0; i<Nn; i++){
		double x = x_min + (x_max - x_min) * i / (Nn - 1);
		for(size_t j = 0; j<Nn; j++){
			double y = y_min + (y_max - y_min) * j / (Nn - 1);
			Br(i, j) =  - evaluate_derivative_Chebyshev_T_expansion(n, Variable::y, a, x, y, x_min, x_max, y_min, y_max) / x;
		}
	}

	dump(argv[2], Br, false);

	// // Compare x derivative
	// for(size_t i = 0; i<N; i++){
	// 	double x = x_min + (x_max - x_min) * i / (N - 1);
	// 	for(size_t j = 0; j<N; j++){
	// 		double y = y_min + (y_max - y_min) * j / (N - 1);
	// 		double exact = df_dx(x, y);
	// 		double cheby = evaluate_derivative_Chebyshev_T_expansion(n, Variable::x, a, x, y, x_min, x_max, y_min, y_max);
	// 		if(std::abs(exact - cheby) > epsilon){
	// 			std::cerr << "d(3 * sin(x) * cos(2 * y + 1))/dx via Chebyshev expansion at (" << x << ", " << y << ") gave " << cheby << " but " << exact << " was expected.\n";
	// 			return 1;
	// 		}
	// 	}
	// }

	// // Compare y derivative
	// for(size_t i = 0; i<N; i++){
	// 	double x = x_min + (x_max - x_min) * i / (N - 1);
	// 	for(size_t j = 0; j<N; j++){
	// 		double y = y_min + (y_max - y_min) * j / (N - 1);
	// 		double exact = df_dy(x, y);
	// 		double cheby = evaluate_derivative_Chebyshev_T_expansion(n, Variable::y, a, x, y, x_min, x_max, y_min, y_max);
	// 		if(std::abs(exact - cheby) > epsilon){
	// 			std::cerr << "d(3 * sin(x) * cos(2 * y + 1))/dy via Chebyshev expansion at (" << x << ", " << y << ") gave " << cheby << " but " << exact << " was expected.\n";
	// 			return 1;
	// 		}
	// 	}
	// }

	return 0;
}

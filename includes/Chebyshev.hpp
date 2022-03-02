/**
 * Chebyshev polynomials are used to express functions
 * in the range [-1, 1]
 * (https://en.wikipedia.org/wiki/Chebyshev_polynomials).
 */
#if !defined(FOCUS_CHEBYSHEV_HPP)
#define FOCUS_CHEBYSHEV_HPP

#include <iostream>
#include <cmath>

#include "types.hpp"
#include "util.hpp"
#include "interpolations.hpp"

/**
 * Calculate Chebyshev polynomials of first kind of
 * order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 */
double Chebyshev_T(int n, double x){
	if (x > 1 || x < -1) {
		std::cerr << "Attempted evaluating T_n(x) out of the domain [-1, 1]\n";
		return nan("");
	}
	return cos(n * acos(x));
}

/**
 * Calculate an array fo Chebyshev polynomials of first kind 
 * from order 0 to order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 * @param T array with n spaces to be filled with the polynomials.
 */
void Chebyshev_T(int n, double x, double T[]){
	if (x > 1 || x < -1) {
		std::cerr << "Attempted evaluating T_n(x) out of the domain [-1, 1]\n";
		return;
	}
	if (n < 0) {
		std::cerr << "Attempted to find Chebyshev polynomials in range but a negative value of n was given.\n";
		return;
	}
	// base cases
	T[0] = 1; if(n > 0) T[1] = x;
	// recurrence relation
	for(int i = 2; i <= n; i++)
		T[i] = 2 * x * T[i - 1] - T[i -2];
}

/**
 * Calculate the singles coefficient of Chebyshev T 
 * expansion for a two variables function from a 
 * matrix of its values.
 * @param M the matrix of values.
 * @param idx first index of coefficient.
 * @param idy second index of coefficient.
 * @param Nx matrix dimension in first variable.
 * @param Ny matrix dimension in second variable.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the coefficient.
 */
double Chebyshev_T_expansion_coefficient(const Matrix2D& M, int idx, int idy, int Nx, int Ny, double x_min, double x_max, double y_min, double y_max){
	double a = 0;
	for (int i = 0; i < Nx; i++){
		double xi = cos((i + 0.5) * PI / Nx);
		double x = x_min + 0.5 * xi * (x_max - x_min);
		for (int j = 0; j < Ny; j++){
			double yj = cos((j + 0.5) * PI / Ny);
			double y = y_min + 0.5 * yj * (y_max - y_min);
			
			a += six_point_formula(x, y, M, x_min, x_max, y_min, y_max) * Chebyshev_T(idx, xi) * Chebyshev_T(idy, yj);
		}
	}
	return a;
}

#if defined(COMPILE_CUDA)

// Add here the kernel to compute all the coefficients in parallell

#else // if defined(COMPILE_CUDA)

/**
 * Calculate the coefficients of Chebyshev T expansion
 * for a two variables function from a matrix of its
 * values to order n.
 * @param n order of expansion.
 * @param a matrix to store the coefficients.
 * @param M the matrix of values.
 * @param Nx matrix dimension in first variable.
 * @param Ny matrix dimension in second variable.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the coefficient.
 */
void Chebyshev_T_expansion(int n, Matrix2D& a, const Matrix2D& M, int Nx, int Ny, double x_min, double x_max, double y_min, double y_max){
	for(int idx = 0; idx < Nx; idx++)
		for(int idy = 0; idy < Ny; idy++)
			a[idx][idy] = Chebyshev_T_expansion_coefficient(M, idx, idy, Nx, Ny, x_min, x_max, y_min, y_max);
}

#endif // if !defined(COMPILE_CUDA)

#endif // FOCUS_CHEBYSHEV_HPP

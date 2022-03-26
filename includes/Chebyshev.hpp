/**
 * Chebyshev polynomials are used to express functions
 * in the range [-1, 1]
 * (https://en.wikipedia.org/wiki/Chebyshev_polynomials).
 */
#if !defined(FOCUS_CHEBYSHEV_HPP)
#define FOCUS_CHEBYSHEV_HPP

#include <iostream>
#include <cmath>

#include "types/matrix2D.hpp"
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
 * Calculate an array of Chebyshev polynomials of first kind 
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
		T[i] = 2 * x * T[i - 1] - T[i - 2];
}

/**
 * Calculate an array of the derivatives of Chebyshev polynomials 
 * of first kind from order 0 to order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 * @param dT array with n spaces to be filled with the polynomials.
 */
void derivative_Chebyshev_T(int n, double x, double dT[]){
	double T[n + 2];
	Chebyshev_T(n + 1, x, T);

	dT[0] = 0;
	for(int i = 1; i <= n; i++)
		dT[i] = 0.5 * i * (T[i - 1] - T[i + 1]) / (1 - x * x);
}

/**
 * Calculate the singles coefficient of Chebyshev T 
 * expansion for a two variables function from a 
 * matrix of its values.
 * @param M the matrix of values.
 * @param n order of expansion.
 * @param idx first index of coefficient.
 * @param idy second index of coefficient.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the coefficient.
 */
double Chebyshev_T_expansion_coefficient(const Matrix2D<double>& M, size_t n, size_t idx, size_t idy, double x_min, double x_max, double y_min, double y_max){
	size_t Nx, Ny;
	Nx = M.shape().first;
	Ny = M.shape().second;
	
	double a = 0;
	for (size_t i = 0; i < Nx; i++){
		double xi = cos((i + 0.5) * pi / Nx);
		for (size_t j = 0; j < Ny; j++){
			double yj = cos((j + 0.5) * pi / Ny);
			a += six_point_formula(xi, yj, M, -1, 1, -1, 1) * Chebyshev_T(idx, xi) * Chebyshev_T(idy, yj);
		}
	}

	// Normalization
	a *= (2.0 - (idx == 0 ? 1 : 0)) / Nx;
	a *= (2.0 - (idy == 0 ? 1 : 0)) / Ny;

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
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 */
void Chebyshev_T_expansion(size_t n, Matrix2D<double>& a, const Matrix2D<double>& M, double x_min, double x_max, double y_min, double y_max){		
	for(size_t idx = 0; idx <= n; idx++)
		for(size_t idy = 0; idy <= n; idy++){
			a(idx, idy) = Chebyshev_T_expansion_coefficient(M, n, idx, idy, x_min, x_max, y_min, y_max);
		}
}

#endif // if !defined(COMPILE_CUDA)

/**
 * Evaluate a two variable function from its
 * Chebyshev expansion coefficients
 * @param n order of expansion
 * @param a matrix of coefficients.
 * @param x first variable value.
 * @param y second variable value.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the function evaluated at (x, y)
 */
double evaluate_Chebyshev_T_expansion(size_t n, const Matrix2D<double>& a, double x, double y, double x_min, double x_max, double y_min, double y_max){
	double v = 0;

	// Normalized to range (-1, 1)
	double xi = (2 * x - (x_min + x_max)) / (x_max - x_min); 
	double yi = (2 * y - (y_min + y_max)) / (y_max - y_min); 

	double Tx[n + 1], Ty[n + 1];
	Chebyshev_T(n, xi, Tx);
	Chebyshev_T(n, yi, Ty);

	for(size_t idx = 0; idx <= n; idx++)
		for(size_t idy = 0; idy <= n; idy++)
			v += a(idx, idy) * Tx[idx] * Ty[idy];

	return v;
}

#endif // FOCUS_CHEBYSHEV_HPP

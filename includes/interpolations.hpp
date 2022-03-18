#if !defined(FOCUS_INTERPOLATIONS_HPP)
#define FOCUS_INTERPOLATIONS_HPP

#include <iostream>
#include <cmath>

#include "types/matrix2D.hpp"

/**
 * Interpolate 2D using the 4 points formula.
 * 
 * Abramowitz, M., Stegun, I. A., & Romer, R. H. (1988). 
 * Handbook of mathematical functions with formulas, graphs,
 * and mathematical tables. 25.2.66
 * 
 * @param x from (x, y) desidered interpolation point.
 * @param y from (x, y) desidered interpolation point.
 * @param f matrix with values of f at equispaced points.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return interpolated f(x, y).
 */
double four_point_formula(double x, double y, const Matrix2D<double>& f, double x_min, double x_max, double y_min, double y_max){
	// The interpolation formula is
	// f(x0 + ph, y0 + qk) = ....
	// where (h, k) are the periods and (x0, y0)
	// is a point of the 2D lattice.

	if (x < x_min || x > x_max || y < y_min || y > y_max){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	int Nx, Ny;
	Nx = f.shape().first;
	Ny = f.shape().second;

	double h = (x_max - x_min) / Nx;
	double k = (y_max - y_min) / Ny;

	size_t i = std::floor((x - x_min) / h);
	size_t j = std::floor((y - y_min) / k);

	double p = x - (x_min + i * h);
	double q = y - (y_min + j * k);

	return (1 - p) * (1 - q) * f(i, j) + p * (1 - q) * f(i + 1, j) + q * (1 - p) * f(i, j + 1) + p * q * f(i + 1, j + 1);
}

/**
 * Interpolate 2D using the 6 points formula.
 * 
 * Abramowitz, M., Stegun, I. A., & Romer, R. H. (1988). 
 * Handbook of mathematical functions with formulas, graphs,
 * and mathematical tables. 25.2.67
 * 
 * @param x from (x, y) desidered interpolation point.
 * @param y from (x, y) desidered interpolation point.
 * @param f matrix with values of f at equispaced points.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return interpolated f(x, y).
 */
double six_point_formula(double x, double y, const Matrix2D<double>& f, double x_min, double x_max, double y_min, double y_max){
	// The interpolation formula is
	// f(x0 + ph, y0 + qk) = ....
	// where (h, k) are the periods and (x0, y0)
	// is a point of the 2D lattice.

	if (x <= x_min || x >= x_max || y <= y_min || y >= y_max){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	int Nx, Ny;
	Nx = f.shape().first;
	Ny = f.shape().second;

	double h = (x_max - x_min) / Nx;
	double k = (y_max - y_min) / Ny;

	size_t i = std::floor((x - x_min) / h);
	size_t j = std::floor((y - y_min) / k);

	if (i == 0 || j == 0) 
		return four_point_formula(x, y, f, x_min, x_max, y_min, y_max);

	double p = x - (x_min + i * h);
	double q = y - (y_min + j * k);

	return 0.5 * q * (q - 1) * f(i, j - 1) + 0.5 * p * (p - 1) * f(i - 1, j) + (1 + p * q - p * p - q * q) * f(i, j) + 0.5 * p * (p - 2 * q + 1) * f(i + 1, j) + 0.5 * q * (q - 2 * p + 1) * f(i, j + 1) + p * q * f(i + 1, j + 1);
}

#endif // FOCUS_INTERPOLATIONS_HPP

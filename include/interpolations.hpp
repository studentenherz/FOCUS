#if !defined(FOCUS_INTERPOLATIONS_HPP)
#define FOCUS_INTERPOLATIONS_HPP

#include <iostream>
#include <cmath>

#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"

/**
 * Interpolate 2D using the 4 points formula.
 * 
 * Abramowitz, M., Stegun, I. A., & Romer, R. H. (1988). 
 * Handbook of mathematical functions with formulas, graphs,
 * and mathematical tables. 25.2.66
 * 
 * @param x from (x, y) desidered interpolation point.
 * @param y from (x, y) desidered interpolation point.
 * @param f scalar field given by a matrix and the (x, y) limits of the space it represents
 * @return interpolated f(x, y).
 */
double four_point_formula(double x, double y, ScalarField f){
	// The interpolation formula is
	// f(x0 + ph, y0 + qk) = ....
	// where (h, k) are the periods and (x0, y0)
	// is a point of the 2D lattice.

	if (x < f.x_min || x >= f.x_max || y < f.y_min || y >= f.y_max){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	size_t Nx, Ny;
	Nx = f.M->shape().first;
	Ny = f.M->shape().second;

	double h = (f.x_max - f.x_min) / (Nx - 1);
	double k = (f.y_max - f.y_min) / (Ny - 1);

	size_t i = std::floor((x - f.x_min) / h);
	size_t j = std::floor((y - f.y_min) / k);

	if (i >= Nx - 1 || j >= Ny - 1){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	double p = (x - (f.x_min + i * h)) / h;
	double q = (y - (f.y_min + j * k)) / k;

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
 * @param f scalar field given by a matrix and the (x, y) limits of the space it represents
 * @return interpolated f(x, y).
 */
double six_point_formula(double x, double y, ScalarField f){
	// The interpolation formula is
	// f(x0 + ph, y0 + qk) = ....
	// where (h, k) are the periods and (x0, y0)
	// is a point of the 2D lattice.

	if (x < f.x_min || x >= f.x_max || y < f.y_min || y >= f.y_max){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	size_t Nx, Ny;
	Nx = f.M->shape().first;
	Ny = f.M->shape().second;

	double h = (f.x_max - f.x_min) / (Nx - 1);
	double k = (f.y_max - f.y_min) / (Ny - 1);

	size_t i = std::floor((x - f.x_min) / h);
	size_t j = std::floor((y - f.y_min) / k);

	if (i >= Nx - 1 || j >= Ny - 1){
		// HERE THERE MUST BE A WARNING LOG
		return nan("");
	}

	if (i == 0 || j == 0) 
		return four_point_formula(x, y, f);

	double p = (x - (f.x_min + i * h)) / h;
	double q = (y - (f.y_min + j * k)) / k;

	return 0.5 * q * (q - 1) * f(i, j - 1) + 0.5 * p * (p - 1) * f(i - 1, j) + (1 + p * q - p * p - q * q) * f(i, j) + 0.5 * p * (p - 2 * q + 1) * f(i + 1, j) + 0.5 * q * (q - 2 * p + 1) * f(i, j + 1) + p * q * f(i + 1, j + 1);
}

#endif // FOCUS_INTERPOLATIONS_HPP

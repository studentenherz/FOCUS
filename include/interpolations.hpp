/**
 * @file interpolations.hpp
 * @brief Methods for interpolating experimental data
 */
#if !defined(FOCUS_INCLUDE_INTERPOLATIONS_HPP)
#define FOCUS_INCLUDE_INTERPOLATIONS_HPP

#include <iostream>
#include <cmath>

#include "types/array.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "types/vector.hpp"

/**
 * [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial)
 * with the 3 closest points.
 * @param xi value in wich to interpolate
 * @param f Array of corresponding f(x) values
 * @param x_min x corresponing to first value in f
 * @param x_max x corresponing to last value in f
 * @return f(xi) interpolated from f
 */
double lagrange_interpolation_3(double xi, const Array<double>& f, double x_min, double x_max){
	double d_x = (x_max - x_min) / (f.size() - 1); // Grid period
	size_t index = std::floor((xi - x_min) / d_x); // Index of closest element to xi to the left

	// Keep index inside boundaries
	if (index <= 0) index = 1;
	if (index >= f.size() - 1) index = f.size() - 2;

	Vector3 x, y, l;
	// Get the points
	for(size_t k = 0; k < 3; k++){
		x[k] = x_min + (index + k - 1) * d_x;
		y[k] = f[index + k  - 1];
	}

	// Basis polynomials evaluated at xi
	for(size_t k = 0; k < 3; k++)
		l[k] = (xi - x[(k + 2) % 3])/(x[k] - x[(k + 2) % 3]) * (xi - x[(k + 1) % 3])/(x[k] - x[(k + 1) % 3]);
				
	// Linear combination
	return dot(y, l);
}

/**
 * [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial)
 * with the 3 closest points.
 * @param xi value in wich to interpolate
 * @param xs Array of x values (sorted ascendingly)
 * @param f Array of corresponding f(x) values
 * @return f(xi) interpolated from f
 */
double lagrange_interpolation_3(double xi, const Array<double>& xs, const Array<double>& f){
	// Index of closest element to xi to the left
	size_t index = 1;
	while (xi > xs[index] && index < xs.size() - 2) index++;
	
	Vector3 x, y, l;
	// Get the points
	for(size_t k = 0; k < 3; k++){
		x[k] = xs[index + k - 1];
		y[k] = f[index + k  - 1];
	}

	// Basis polynomials evaluated at xi
	for(size_t k = 0; k < 3; k++)
		l[k] = (xi - x[(k + 2) % 3])/(x[k] - x[(k + 2) % 3]) * (xi - x[(k + 1) % 3])/(x[k] - x[(k + 1) % 3]);
				
	// Linear combination
	return dot(y, l);
}


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
	Nx = f.M.shape().first;
	Ny = f.M.shape().second;

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
	Nx = f.M.shape().first;
	Ny = f.M.shape().second;

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

#endif // FOCUS_INCLUDE_INTERPOLATIONS_HPP

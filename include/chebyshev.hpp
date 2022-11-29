/**
 * @file chebyshev.hpp
 * @brief [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
 * are used in this program in order to express functions in the range [-1, 1].
 * 
 * This part of the code deals with part of calculating the magnetic field from the 
 * poloidal flux \f$\psi(r, z)\f$ and current \f$ I(r, z)\f$
 * \f[
 * 	\vec{B} = \frac{1}{r}\left(\frac{\partial \psi}{\partial r} \hat{z} -\frac{\partial \psi}{\partial z} \hat{r} \right) + \frac{I}{r} \hat{\theta}
 * \f]
 * 
 * The magnetic field needs to have divergence zero \f$\nabla \cdot \vec{B} = 0 \f$ 
 * in order to be compliant with Maxwell's equations, but the numerical nature of the data attempts 
 * against it. In that regard, an expansion of the poloidal flux in the form \f$\psi(r, z) = R(r)Z(z)\f$
 * will ensure the null divergence by construction (keeping in mind that \f$ d/d\theta \equiv 0 \f$ because 
 * of the ax-symmetry). 
 * 
 * @ref ChebyshevExpansion implements such an expansion using first kind Chebyshev polynomials 
 * \f[
 * 	\psi(r, z) = \sum_{i, j}^{n}a_{i, j}T_i(x(r))T_j(y(z))
 * \f]
 * where
 * \f[
 * 	a_{i, j} = \left(\frac{2 - \delta_{0k}}{N_x}\right) \left(\frac{2 - \delta_{0l}}{N_y} \right) \sum_{i}^{N_x}\sum_{j}^{N_y} \Psi(r_k, z_l)T_i(x_k)T_j(y_l)
 * \f]
 * and
 * \f[
 * 	\begin{align}
 *  	x_k &= \cos\left((k + 1/2)\pi/N_x\right)\\
 *  	y_l &= \cos\left((l + 1/2)\pi/N_y\right) 
 *	\end{align}
 * \f]
 * where \f$ N_x\f$ and \f$ N_y\f$ are the dimensions of the matrix of experimental data.
 */
#if !defined(FOCUS_INCLUDE_CHEBYSHEV_HPP)
#define FOCUS_INCLUDE_CHEBYSHEV_HPP

#include <iostream>
#include <cmath>

#include "types/array.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "util.hpp"
#include "interpolations.hpp"

/**
 * Calculate Chebyshev polynomials of first kind of
 * order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 */
double Chebyshev_T(size_t n, double x){
	if (x > 1 || x < -1) {
		std::cerr << "Attempted evaluating T_n(" << x << ") out of the domain [-1, 1]\n";
		return nan("");
	}
	return cos(n * acos(x));
}

/**
 * Calculate Chebyshev polynomials of second kind of
 * order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 */
double Chebyshev_U(size_t n, double x){
	if (x > 1 || x < -1) {
		std::cerr << "Attempted evaluating U_n(" << x << ") out of the domain [-1, 1]\n";
		return nan("");
	}
	return sin((n + 1) * acos(x)) / sin(acos(x));
}

/**
 * Calculate an array of Chebyshev polynomials of first kind 
 * from order 0 to order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 * @param T array with n spaces to be filled with the polynomials.
 * @return true if no error occurred.
 */
bool Chebyshev_T(size_t n, double x, Array<double>& T){
	if(x > 1 || x < -1){
		std::cerr << "Attempted evaluating T_n(" << x << ") out of the domain [-1, 1]\n";		
		for (size_t i = 0; i <= n; i++)
			T[i] = nan("");
		return false;
	}
	// base cases
	T[0] = 1; if(n > 0) T[1] = x;
	// recurrence relation
	for(size_t i = 2; i <= n; i++)
		T[i] = 2 * x * T[i - 1] - T[i - 2];
	
	return true;
}

/**
 * Calculate an array of Chebyshev polynomials of second kind 
 * from order 0 to order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 * @param U array with n spaces to be filled with the polynomials.
 * @return true if no error occurred.
 */
bool Chebyshev_U(size_t n, double x, Array<double>& U){
	if(x > 1 || x < -1){
		std::cerr << "Attempted evaluating T_n(" << x << ") out of the domain [-1, 1]\n";		
		for (size_t i = 0; i <= n; i++)
			U[i] = nan("");
		return false;
	}
	// base cases
	U[0] = 1; if(n > 0) U[1] = 2 * x;
	// recurrence relation
	for(size_t i = 2; i <= n; i++)
		U[i] = 2 * x * U[i - 1] - U[i - 2];
	
	return true;
}

/**
 * Calculate an array of the derivatives of Chebyshev polynomials 
 * of first kind from order 0 to order n at x.
 * @param n order to which calculate the polynomials.
 * @param x value to evaluate the polynomials.
 * @param dT array with n spaces to be filled with the polynomials.
 * @return true if no errors occurred.
 */
bool derivative_Chebyshev_T(size_t n, double x, Array<double>& dT){
	// https://en.wikipedia.org/wiki/Chebyshev_polynomials#Differentiation_and_integration
	Array<double> U(n);
	if (!Chebyshev_U(n - 1, x, U)){
		for(size_t i = 0; i <= n; i++)
			dT[i] = nan("");

		return false;
	}

	dT[0] = 0;
	for(size_t i = 1; i <= n; i++)
		dT[i] = i * U[i - 1];

	return true;
}

/**
 * Calculate the singles coefficient of Chebyshev T 
 * expansion for a two variables function from a 
 * matrix of its values.
 * @param M the matrix of values.
 * @param idx first index of coefficient.
 * @param idy second index of coefficient.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the coefficient.
 */
double Chebyshev_T_expansion_coefficient(size_t idx, size_t idy, ScalarField f, double x_min, double x_max, double y_min, double y_max){
	size_t Nx, Ny;
	Nx = f.M.shape().first;
	Ny = f.M.shape().second;
	
	double x_m = 0.5 * (x_min + x_max);
	double y_m = 0.5 * (y_min + y_max);
	double Dx = 0.5 * (x_max - x_min);
	double Dy = 0.5 * (y_max - y_min);

	double a = 0;
	for (size_t i = 0; i < Nx; i++){
		double xi = cos((i + 0.5) * pi / Nx);
		double x = x_m + xi * Dx;
		for (size_t j = 0; j < Ny; j++){
			double yj = cos((j + 0.5) * pi / Ny);
			double y = y_m + yj * Dy;
			a += six_point_formula(x, y, f) * Chebyshev_T(idx, xi) * Chebyshev_T(idy, yj);
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
void Chebyshev_T_expansion(size_t n, Matrix2D<double>& a, ScalarField f, double x_min, double x_max, double y_min, double y_max){		
	for(size_t idx = 0; idx <= n; idx++)
		for(size_t idy = 0; idy <= n; idy++){
			a(idx, idy) = Chebyshev_T_expansion_coefficient(idx, idy, f, x_min, x_max, y_min, y_max);
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

	Array<double> Tx(n + 1), Ty(n + 1);
	Chebyshev_T(n, xi, Tx);
	Chebyshev_T(n, yi, Ty);

	for(size_t idx = 0; idx <= n; idx++)
		for(size_t idy = 0; idy <= n; idy++)
			v += a(idx, idy) * Tx[idx] * Ty[idy];

	return v;
}

/** Enumeration of possible variables
 */
enum Variable {x, y};

/**
 * Evaluate a two variable function's derivative 
 * from its Chebyshev expansion coefficients.
 * @param n order of expansion
 * @param var variable with respect to which calculate the derivative
 * @param a matrix of coefficients.
 * @param x first variable value.
 * @param y second variable value.
 * @param x_min minimum value of x represented in the matrix.
 * @param x_max maximum value of x represented in the matrix.
 * @param y_min minimum value of y represented in the matrix.
 * @param y_max maximum value of y represented in the matrix.
 * @return the function evaluated at (x, y)
 */
double evaluate_derivative_Chebyshev_T_expansion(size_t n, Variable var, const Matrix2D<double>& a, double x, double y, double x_min, double x_max, double y_min, double y_max){
	double v = 0;

	// Normalized to range (-1, 1)
	double xi = (2 * x - (x_min + x_max)) / (x_max - x_min); 
	double yi = (2 * y - (y_min + y_max)) / (y_max - y_min); 

	Array<double> T(n + 1), dT(n + 1);

	if(var == Variable::x){
		derivative_Chebyshev_T(n, xi, dT);
		Chebyshev_T(n, yi, T);
	}
	else{
		Chebyshev_T(n, xi, T);
		derivative_Chebyshev_T(n, yi, dT);
	}

	for(size_t idx = 0; idx <= n; idx++)
		for(size_t idy = 0; idy <= n; idy++)
			v += a(idx, idy) * (var == Variable::x ? dT[idx] * T[idy] : T[idx] * dT[idy]);

	// Normalization
	v *= 2 / (var == Variable::x ? (x_max - x_min) : (y_max - y_min));

	return v;
}

/**
 * A class to encapsulate all Chebyshev expansion related functions
 */ 
class ChebyshevExpansion{
	size_t n;	// order of expansion
	Matrix2D<double> a; // matrix of coefficients
	double x_min, x_max, y_min, y_max; // limits
public:

	/**
	 * Calculates the Chebyshev expansion of the scalar field
	 * up to order n in the given x, y region
	 * @param n order of the expansion
	 * @param f scalar field to be expanded
	 * @param x_min minimum value of x of the expansion
	 * @param x_max maximum value of x of the expansion
	 * @param y_min minimum value of y of the expansion
 	 * @param y_max maximum value of y of the expansion
	 */
	ChebyshevExpansion(size_t order, ScalarField f, double xmin, double xmax, double ymin, double ymax): n(order), a(n + 1, n + 1) {
		/* Adding extra room for x and y so that Chebyshev is not
			evaluated in 1 or -1; depending on the numbers, the floating
			representation might give something like 1.00000000005 that 
			still breaks the function.
		*/

		double ex = 0.01 * (xmax - xmin);
		double ey = 0.01 * (ymax - ymin);

		x_min = std::max(f.x_min, xmin - ex);
		x_max = std::min(f.x_max, xmax + ex);
		y_min = std::max(f.y_min, ymin - ey);
		y_max = std::min(f.y_max, ymax + ey);

		// Calculate the matrix of coefficients
		Chebyshev_T_expansion(n, a, f, x_min, x_max, y_min, y_max);
	}

	/**
	 * Evaluate the scalar field f from the expansion coefficients
	 * @param x x
	 * @param y y
	 * @return f(x, y)
	 */
	double operator()(double x, double y){
		return evaluate_Chebyshev_T_expansion(n, a, x, y, x_min, x_max, y_min, y_max);
	}

	/**
	 * Evaluate the derivative with respect to x of the 
	 * scalar field f from the expansion coefficients
	 * @param x x
	 * @param y y
	 * @return df_dx(x, y)
	 */
	double dx(double x, double y){
		return evaluate_derivative_Chebyshev_T_expansion(n, Variable::x, a, x, y, x_min, x_max, y_min, y_max);
	}

	/**
	 * Evaluate the derivative with respect to y of the 
	 * scalar field f from the expansion coefficients
	 * @param x x
	 * @param y y
	 * @return df_dy(x, y)
	 */
	double dy(double x, double y){
		return evaluate_derivative_Chebyshev_T_expansion(n, Variable::y, a, x, y, x_min, x_max, y_min, y_max);
	}
};

#endif // FOCUS_INCLUDE_CHEBYSHEV_HPP

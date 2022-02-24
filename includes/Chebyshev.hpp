/**
 * Chebyshev polynomials are used to express functions
 * in the range [-1, 1]
 * (https://en.wikipedia.org/wiki/Chebyshev_polynomials).
 */
#if !defined(FOCUS_CHEBYSHEV_HPP)
#define FOCUS_CHEBYSHEV_HPP

#include <iostream>
#include <cmath>

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
	for(int i = 2; i <= n; i++)
		T[i] = 2 * x * T[i - 1] - T[i -2];
}

#endif // FOCUS_CHEBYSHEV_HPP

#if !defined(FOCUS_INCLUDE_TYPES_PLASMA_HPP)
#define FOCUS_INCLUDE_TYPES_PLASMA_HPP

#include "types/array.hpp"
#include "types/matrix_2d.hpp"

struct Plasma{
	int shot;
	size_t nexp;
	size_t nion;
	double masse;
	double ze;

	Array<double> mass;
	Array<double> z;

	Array<double> psi;

	Array<double> ne;
	Matrix2D<double> ni;	

	Array<double> te;
	Matrix2D<double> ti;	

	Plasma(int shot, size_t nexp, size_t nion): shot(shot), nexp(nexp), nion(nion), mass(nion), z(nion), psi(nexp), ne(nexp), ni(nion, nexp), te(nexp), ti(nion, nexp) {} 
};

#endif // FOCUS_INCLUDE_TYPES_PLASMA_HPP

#if !defined(FOCUS_INCLUDE_TYPES_PLASMA_HPP)
#define FOCUS_INCLUDE_TYPES_PLASMA_HPP

#include "types/array.hpp"

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
	Array<Array<double>> ni;	

	Array<double> te;
	Array<Array<double>> ti;	

	Plasma(int shot, size_t nexp, size_t nion): shot(shot), nexp(nexp), nion(nion), mass(nion), z(nion), psi(nexp), ne(nexp), ni(nion), te(nexp), ti(nion) {
		for (size_t i = 0; i < nion; i++){
			ni[i].resize(nexp);
			ti[i].resize(nexp);
		}
	} 
};

#endif // FOCUS_INCLUDE_TYPES_PLASMA_HPP

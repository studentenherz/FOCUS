#if !defined(FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP)
#define FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP

#include "types/matrix_2d.hpp"

struct ScalarField{
	Matrix2D<double>& M;
	double x_min;
	double x_max;
	double y_min;
	double y_max;

	ScalarField(Matrix2D<double>& _M, double _x_min, double _x_max, double _y_min, double _y_max): M(_M), x_min(_x_min), x_max(_x_max), y_min(_y_min), y_max(_y_max) {}

	double operator()(int i, int j) const{
		return M(i, j);
	}
};

#endif // FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP



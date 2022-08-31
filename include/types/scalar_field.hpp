#if !defined(FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP)
#define FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP

#include "types/matrix_2d.hpp"

/**
 * Wrapper that encapsulates a matrix description of
 * a scalar field together with the limits of the space
 * it represents
 */
struct ScalarField{
	Matrix2D<double>& M;	///< Matrix description of the field

	/**
	 * @name Space limits
	 */
	///@{
	double x_min;		///< minimum x
	double x_max;		///< maximum x
	double y_min;		///< minimum y
	double y_max;		///< maximum y
	///@}

	/**
	 * Constructor
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	ScalarField(Matrix2D<double>& _M, double _x_min, double _x_max, double _y_min, double _y_max): M(_M), x_min(_x_min), x_max(_x_max), y_min(_y_min), y_max(_y_max) {}

	/**
	 * Propagate the matrix access operator
	 * @param i first index
	 * @param j second index
	 * @return `M(i, j)`
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	double operator()(int i, int j) const {
		return M(i, j);
	}
};

#endif // FOCUS_INCLUDE_TYPES_SCALAR_FIELD_HPP



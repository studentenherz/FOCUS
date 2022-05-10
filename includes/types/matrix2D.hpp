#if !defined(FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP)
#define FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP

#include "pair.hpp"

using shape_t = pair<size_t>;

template <typename T>
class Matrix2D{
	T *arr;
	shape_t _shape;
public:
	/**
	 * Default constructor
	 * @param n size on the first axis
	 * @param m size on the second axis
	 */
	Matrix2D(size_t n = 0, size_t m = 0){
		_shape = {n, m};
		arr = new T[n * m + 1];
	}

	/**
	 * Move constructor; this allows the matrix to be passed
	 * as a return value from a function sapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	Matrix2D(Matrix2D&& other): _shape(other._shape){
		arr = other.arr;
		other.arr = new T[1];
	}

	/**
	 * Desctructor; dealocate heap and set pointer to null
	 */
	~Matrix2D(){
		delete[] arr;
		arr = NULL;
	}

	/** 
	 * Change shape of matrix droping previous data
	 * @param n size on the first axis
	 * @param m size on the second axis
	 */ 
	void reshape(size_t n, size_t m){
		delete[] arr;
		arr = new T[n * m + 1];
		_shape = {n, m};
	}

	/**
	 * Get shape of the matrix
	 * @return shape
	 */
	shape_t shape() const {
		return _shape;
	}

	/** 
	 * Write access to matrix
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	T &operator()(size_t i, size_t j){
		if (i > _shape.first || j > _shape.second)
			return arr[_shape.first * _shape.second];
		return arr[i * _shape.second + j];
	}
	
	/** 
	 * Read-only access to matrix
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	const T &operator()(size_t i, size_t j) const {
		if (i > _shape.first || j > _shape.second)
			return arr[_shape.first * _shape.second];
		return arr[i * _shape.second + j];
	}
};

#endif // FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP

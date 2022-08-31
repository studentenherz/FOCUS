#if !defined(FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP)
#define FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP

#include "pair.hpp"

using shape_t = pair<size_t>;

template <typename T>
class Matrix2D{
	T *_arr;
	shape_t _shape;
	bool _copied;
public:
	/**
	 * Default constructor
	 * @param n size on the first axis
	 * @param m size on the second axis
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Matrix2D(size_t n = 0, size_t m = 0){
		_copied = false;
		_shape = {n, m};
		_arr = new T[n * m + 1];
	}

	Matrix2D& operator=(const Matrix2D&&) = delete;

	/**
	 * Copy constructor; put copied flag
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Matrix2D(Matrix2D& other){
		_shape = other._shape;
		_arr = other._arr;
		_copied = true;
	}

	/**
	 * Move constructor; this allows the matrix to be passed
	 * as a return value from a function swapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Matrix2D(Matrix2D&& other){
		_shape = other._shape;
		_arr = other._arr;
		_copied = other._copied;

		// If the previous matrix was a copied one just copy
		// else actually do the move
		if (!_copied)
			other._arr = NULL;
	}

	/**
	 * Construct in host for device from Matrix2D 
	 * @param other matrix2D to construct from
	 */
	#ifdef __CUDACC__
	__host__
	void construct_in_host_for_device(Matrix2D<T>& other){
		_copied = true;
		_shape = other._shape;
		size_t _size = _shape.first * _shape.second;
		cudaMalloc(&_arr, sizeof(T) * (_size + 1));
		cudaMemcpy(_arr, other._arr, sizeof(T) * _size, cudaMemcpyHostToDevice);
	}
	#endif

	/**
	 * Desctructor; dealocate heap and set pointer to null
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	~Matrix2D(){
		if (_copied) return;
		delete[] _arr;
		_arr = NULL;
	}

	/** 
	 * Change shape of matrix droping previous data
	 * @param n size on the first axis
	 * @param m size on the second axis
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif 
	void reshape(size_t n, size_t m){
		delete[] _arr;
		_arr = new T[n * m + 1];
		_shape = {n, m};
	}

	/**
	 * Get shape of the matrix
	 * @return shape
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	shape_t shape() const {
		return _shape;
	}

	/** 
	 * Write access to matrix
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	T &operator()(size_t i, size_t j){
		if (i > _shape.first || j > _shape.second)
			return _arr[_shape.first * _shape.second];
		return _arr[i * _shape.second + j];
	}
	
	/** 
	 * Read-only access to matrix
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	const T &operator()(size_t i, size_t j) const {
		if (i > _shape.first || j > _shape.second)
			return _arr[_shape.first * _shape.second];
		return _arr[i * _shape.second + j];
	}
};

#endif // FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP

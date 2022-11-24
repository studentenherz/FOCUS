#if !defined(FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP)
#define FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP

#include "pair.hpp"
#include "handle_cuda_errors.hpp"

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
		_arr = new T[n * m];
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
	cudaError_t construct_in_host_for_device(Matrix2D<T>& other){
		_copied = true;
		_shape = other._shape;
		size_t _size = _shape.first * _shape.second;
		propagateCUDAErr( cudaMalloc(&_arr, sizeof(T) * _size) );
		propagateCUDAErr( cudaMemcpy(_arr, other._arr, sizeof(T) * _size, cudaMemcpyHostToDevice) );
		return cudaSuccess;
	}
	#endif

	/**
	 * Copy to host from device
	 * @param other matrix2D to construct from
	 */
	#ifdef __CUDACC__
	__host__
	cudaError_t copy_to_host_from_device(Matrix2D<T>& other){
		_copied = false;
		_shape = other._shape;
		size_t _size = _shape.first * _shape.second;
		propagateCUDAErr( cudaMemcpy(_arr, other._arr, sizeof(T) * _size, cudaMemcpyDeviceToHost) );
		return cudaSuccess;
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
		_arr = new T[n * m];
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

	#ifdef __CUDA_ARCH__

	/** 
	 * Write access to matrix (device)
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	__device__
	T &operator()(size_t i, size_t j){
		if (i >= _shape.first || j >= _shape.second){
			printf("Error, out of bounds access to Matrix2D %s %d\n", __FILE__, __LINE__);
			assert(0);
		}
		return _arr[i * _shape.second + j];
	}
	
	/** 
	 * Read-only access to matrix (device)
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	__device__
	const T &operator()(size_t i, size_t j) const {
		if (i >= _shape.first || j >= _shape.second){
			printf("Error, out of bounds access to Matrix2D %s %d\n", __FILE__, __LINE__);
			assert(0);
		}
		return _arr[i * _shape.second + j];
	}

	#else

	/** 
	 * Write access to matrix (host)
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	#ifdef __CUDACC__
	__host__
	#endif
	T &operator()(size_t i, size_t j){
		if (i >= _shape.first || j >= _shape.second){
			fprintf(stderr, "Error, out of bounds access to Matrix2D %s %d\n", __FILE__, __LINE__);
			exit(2);
		}
		return _arr[i * _shape.second + j];
	}
	
	/** 
	 * Read-only access to matrix (host)
	 * @param i index on the first axis
	 * @param j index on the second axis
	 * @return M(i, j)
	 */
	#ifdef __CUDACC__
	__host__
	#endif
	const T &operator()(size_t i, size_t j) const {
		if (i >= _shape.first || j >= _shape.second){
			fprintf(stderr, "Error, out of bounds access to Matrix2D %s %d\n", __FILE__, __LINE__);
			exit(2);
		}
		return _arr[i * _shape.second + j];
	}

	#endif // __CUDA_ARCH__
};

#endif // FOCUS_INCLUDE_TYPES_MATRIX_2D_HPP

#if !defined(FOCUS_INCLUDE_TYPES_ARRAY_HPP)
#define FOCUS_INCLUDE_TYPES_ARRAY_HPP

#include "handle_cuda_errors.hpp"

template <typename T>
class Array{
	T *_arr;
	size_t _size;
	bool _copied;
public:
	/**
	 * Default constructor
	 * @param n size of array
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Array(size_t n = 0): _size(n), _copied(false) {
		_arr = new T[_size];
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Array(std::initializer_list<T> l): _size(l.size()) {
		_arr = new T[_size];
    _copied = false;

		size_t index = 0;
		for(const T* it = l.begin(); it != l.end(); it++)
			_arr[index++] = *it;
	}

	/**
	 * Construct in host for device from array
	 * @param other_arr array to construct from
	 * @param n size of array
	 */
	#ifdef __CUDACC__
	__host__
	Array(T* other_arr, size_t size): _size(size), _copied(false) {
		gpuErrchk(cudaMalloc(&_arr, sizeof(T) * (_size)));
		gpuErrchk(cudaMemcpy(_arr, other_arr, sizeof(T) * _size, cudaMemcpyHostToDevice));
	}
	#endif

	/**
	 * Move constructor; this allows the array to be passed
	 * as a return value from a function sapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Array(Array&& other){
		_size = other._size;
		_arr = other._arr;
		_copied = other._copied;

		// If the previous array was a copied one just copy
		// else actually do the move
		if (!_copied)
			other._arr = NULL;
	}

	/**
	 * Copy constructor; declaring the move constructor makes
	 * compiler implicitly declare the copy constructor deleted
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Array(Array& other){
		_copied = true;
		_size = other._size;
		_arr = other._arr;
	}

	/**
	 * Construct in host for device from Array 
	 * @param other Array to construct from
	 */
	#ifdef __CUDACC__
	__host__
	cudaError_t construct_in_host_for_device(Array<T>& other){
		_copied = true;
		_size = other._size;
		propagateCUDAErr( cudaMalloc(&_arr, sizeof(T) * _size) );
		propagateCUDAErr( cudaMemcpy(_arr, other._arr, sizeof(T) * _size, cudaMemcpyHostToDevice) );
		return cudaSuccess;
	}
	#endif
	
	/**
	 * Copy to host from device
	 * @param other Array to construct from
	 */
	#ifdef __CUDACC__
	__host__
	cudaError_t copy_to_host_from_device(Array<T>& other){
		_copied = false;
		_size = other._size;
		_arr = new T[_size];
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
	~Array(){
		if (_copied) return;
		delete[] _arr;
		_arr = NULL;
	}

	#ifdef __CUDACC__
	__host__
	cudaError_t deallocate_cuda(){
		propagateCUDAErr( cudaFree(_arr) );
		return cudaSuccess;
	}
	#endif

	#ifdef __CUDACC__
	__host__
	T from_device_at(size_t i){
		T h_val;
		gpuErrchk( cudaMemcpy(&h_val, &_arr[i], sizeof(T), cudaMemcpyDeviceToHost) );
		return h_val;
	}
	#endif

	#ifdef __CUDA_ARCH__

	/**
	 * Write access to the array (device)
	 * @param i index
	 * @return reference to i-th element
	 */
	__device__
	T &operator[](size_t i){
		if (i >= _size){
			printf("Error, out of bounds access to Array %s %d\n", __FILE__, __LINE__);
			assert(0);
		}
		return _arr[i];
	}
	
	/**
	 * Read only access to the array (device)
	 * @param i index
	 * @return reference to i-th element
	 */
	__device__
	const T &operator[](size_t i) const {
		if (i >= _size){
			printf("Error, out of bounds access to Array %s %d\n", __FILE__, __LINE__);
			assert(0);
		}
		return _arr[i];
	}

	#else

	/**
	 * Write access to the array (host)
	 * @param i index
	 * @return reference to i-th element
	 */
	#ifdef __CUDACC__
	__host__
	#endif
	T &operator[](size_t i){
		if (i >= _size){
			fprintf(stderr, "Error, out of bounds access to Array %s %d\n", __FILE__, __LINE__);
			exit(2);
		}
		return _arr[i];
	}
	
	/**
	 * Read only access to the array (host)
	 * @param i index
	 * @return reference to i-th element
	 */
	#ifdef __CUDACC__
	__host__
	#endif
	const T &operator[](size_t i) const {
		if (i >= _size){
			fprintf(stderr, "Error, out of bounds access to Array %s %d\n", __FILE__, __LINE__);
			exit(2);
		}
		return _arr[i];
	}

	#endif	// __CUDA_ARCH__

	/** 
	 * Get array size
	 * @return array size
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	size_t size() const {
		return _size;
	}

	/** 
	 * Resize array droping stored values
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void resize(size_t n){
		delete[] _arr;
		_size = n;
		_arr = new T[_size];
	}
}; // class Array

/**
 * Returns the smallest element from an array
 * @param a Array
 * @return smallest element of `a`
 */
template<typename T>
#ifdef __CUDACC__
__host__ __device__
#endif
T min(const Array<T>& a){
	T m = a[0];
	for(size_t i = 1; i < a.size(); i++)
		if (a[i] < m)
			m = a[i];
	return m;
}

/**
 * Returns the larges element from an array
 * @param a Array
 * @return larges element of `a`
 */
template<typename T>
#ifdef __CUDACC__
__host__ __device__
#endif
T max(const Array<T>& a){
	T m = a[0];
	for(size_t i = 1; i < a.size(); i++)
		if (a[i] > m)
			m = a[i];
	return m;
}

#endif // FOCUS_INCLUDE_TYPES_ARRAY_HPP

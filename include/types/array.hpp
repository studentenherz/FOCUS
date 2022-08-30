#if !defined(FOCUS_INCLUDE_TYPES_ARRAY_HPP)
#define FOCUS_INCLUDE_TYPES_ARRAY_HPP

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
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	Array(size_t n = 0): _size(n), _copied(false) {
		_arr = new T[_size + 1];
	}

	/**
	 * Construct in host for device from array
	 * @param other_arr array to construct from
	 * @param n size of array
	 */
	#ifdef CUDA_BUILD
	__host__
	Array(T* other_arr, size_t size, bool from_host = true): _size(size), _copied(from_host) {
		cudaMalloc(&_arr, sizeof(T) * (_size + 1));
		cudaMemcpy(_arr, other_arr, sizeof(T) * _size, cudaMemcpyHostToDevice);
	}
	#endif

	/**
	 * Move constructor; this allows the array to be passed
	 * as a return value from a function sapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	#ifdef CUDA_BUILD
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
	#ifdef CUDA_BUILD
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
	#ifdef CUDA_BUILD
	__host__
	void construct_in_host_for_device_from_array(Array<T>& other){
		_copied = true;
		_size = other._size;
		cudaMalloc(&_arr, sizeof(T) * (_size + 1));
		cudaMemcpy(_arr, other._arr, sizeof(T) * _size, cudaMemcpyHostToDevice);
	}
	#endif

	/**
	 * Desctructor; dealocate heap and set pointer to null
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	~Array(){
		if (_copied) return;
		delete[] _arr;
		_arr = NULL;
	}

	/**
	 * Write access to the array
	 * @param i index
	 * @return reference to i-th element
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	T &operator[](size_t i){
		if (i > _size)
			return _arr[_size];
		return _arr[i];
	}
	
	/**
	 * Read only access to the array
	 * @param i index
	 * @return reference to i-th element
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	const T &operator[](size_t i) const {
		if (i > _size)
			return _arr[_size];
		return _arr[i];
	}

	/** 
	 * Get array size
	 * @return array size
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	size_t size() const {
		return _size;
	}

	/** 
	 * Resize array droping stored values
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	void resize(size_t n){
		delete[] _arr;
		_size = n;
		_arr = new T[_size + 1];
	}
}; // class Array

/**
 * Returns the smallest element from an array
 * @param a Array
 * @return smallest element of `a`
 */
template<typename T>
#ifdef CUDA_BUILD
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
#ifdef CUDA_BUILD
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

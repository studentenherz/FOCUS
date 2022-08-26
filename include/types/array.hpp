#if !defined(FOCUS_INCLUDE_TYPES_ARRAY_HPP)
#define FOCUS_INCLUDE_TYPES_ARRAY_HPP

template <typename T>
class Array{
	T *arr;
	size_t _size;
	bool _in_device_from_host;
public:
	/**
	 * Default constructor
	 * @param n size of array
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	Array(size_t n = 0): _size(n), _in_device_from_host(false) {
		arr = new T[_size + 1];
	}

	/**
	 * Construct from array
	 * @param other_arr array to construct from
	 * @param n size of array
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	Array(T* other_arr, size_t n, bool in_device_from_host = true): _size(n), _in_device_from_host(in_device_from_host) {
		arr = other_arr;
	}
	
	/**
	 * Move constructor; this allows the array to be passed
	 * as a return value from a function sapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	Array(Array&& other): _size(other._size), _in_device_from_host(other._in_device_from_host) {
		arr = other.arr;
		other.arr = NULL;
	}

	/**
	 * Copy constructor; declaring the move constructor makes
	 * compiler implicitly declare the copy constructor deleted
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	Array(Array& other): Array(other._size) {
		for(size_t i = 0; i < _size; i++)
			arr[i] = other.arr[i];
	}

	/**
	 * Desctructor; dealocate heap and set pointer to null
	 */
	#ifdef CUDA_BUILD
	__host__ __device__
	#endif
	~Array(){
		if (_in_device_from_host) return;
		delete[] arr;
		arr = NULL;
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
			return arr[_size];
		return arr[i];
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
			return arr[_size];
		return arr[i];
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
		delete[] arr;
		_size = n;
		arr = new T[_size + 1];
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

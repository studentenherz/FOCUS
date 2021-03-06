#if !defined(FOCUS_INCLUDE_TYPES_ARRAY_HPP)
#define FOCUS_INCLUDE_TYPES_ARRAY_HPP

template <typename T>
class Array{
	T *arr;
	size_t _size;
public:
	/**
	 * Default constructor
	 * @param n size of array
	 */
	Array(size_t n = 0){
		_size = n;
		arr = new T[_size + 1];
	}
	
	/**
	 * Move constructor; this allows the array to be passed
	 * as a return value from a function sapping the pointers
	 * and keeping the allocated data on the heap.
	 */
	Array(Array&& other): _size(other._size) {
		arr = other.arr;
		other.arr = NULL;
	}

	/**
	 * Copy constructor; declaring the move constructor makes
	 * compiler implicitly declare the copy constructor deleted
	 */
	Array(Array& other): Array(other._size) {
		for(size_t i = 0; i < _size; i++)
			arr[i] = other.arr[i];
	}

	/**
	 * Desctructor; dealocate heap and set pointer to null
	 */
	~Array(){
		delete[] arr;
		arr = NULL;
	}

	/**
	 * Write access to the array
	 * @param i index
	 * @return reference to i-th element
	 */
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
	const T &operator[](size_t i) const {
		if (i > _size)
			return arr[_size];
		return arr[i];
	}

	/** 
	 * Get array size
	 * @return array size
	 */
	size_t size() const {
		return _size;
	}

	/** 
	 * Resize array droping stored values
	 */
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
T min(const Array<T>& a){
	T m = a[0];
	for(size_t i = 1; i < a.size(); i++)
		m = std::min(m, a[i]);
	return m;
}

/**
 * Returns the larges element from an array
 * @param a Array
 * @return larges element of `a`
 */
template<typename T>
T max(const Array<T>& a){
	T m = a[0];
	for(size_t i = 1; i < a.size(); i++)
		m = std::max(m, a[i]);
	return m;
}

#endif // FOCUS_INCLUDE_TYPES_ARRAY_HPP

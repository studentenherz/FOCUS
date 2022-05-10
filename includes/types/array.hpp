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
	Array(Array&& other): _size(other._size){
		arr = other.arr;
		other.arr = new T[1];
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
		return arr[i];
	}
	
	/**
	 * Read only access to the array
	 * @param i index
	 * @return reference to i-th element
	 */
	const T &operator[](size_t i) const {
		return arr[i];
	}

	/** 
	 * Get array size
	 * @return array size
	 */
	size_t size(){
		return _size;
	}
};

#endif // FOCUS_INCLUDE_TYPES_ARRAY_HPP

#if !defined(FOCUS_INCLUDE_TYPES_ARRAY_HPP)
#define FOCUS_INCLUDE_TYPES_ARRAY_HPP

template <typename T>
class Array{
	T *arr;
	size_t _size;
public:
	Array(){
		_size = 0;
		arr = new T[1];
	}
	Array(size_t n){
		_size = n;
		arr = new T[_size + 1];
	}
	~Array(){
		delete[] arr;
		arr = NULL;
	}
	T &operator[](size_t i){
		return arr[i];
	}
	const T &operator[](size_t i) const {
		return arr[i];
	}
	size_t size(){
		return _size;
	}
};

#endif // FOCUS_INCLUDE_TYPES_ARRAY_HPP

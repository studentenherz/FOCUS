#if !defined(FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP)
#define FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP

#include "pair.hpp"

using shape_t = pair<size_t>;

template <typename T>
class Matrix2D{
	T *arr;
	shape_t _shape;
public:
	Matrix2D() {
		_shape = {0, 0};
		arr = new T[1];
	}
	Matrix2D(size_t n, size_t m){
		arr = new T[n * m + 1];
		_shape = {n, m};
	}
	Matrix2D(shape_t other_shape) {
		_shape = other_shape;
		arr = new T[_shape.first * _shape.second + 1];
	}
	~Matrix2D(){
		delete[] arr;
		arr = NULL;
	}
	void reshape(size_t n, size_t m){
		delete[] arr;
		arr = new T[n * m + 1];
		_shape = {n, m};
	}
	shape_t shape() const {
		return _shape;
	}
	T &operator()(size_t i, size_t j){
		if (i > _shape.first || j > _shape.second)
			return arr[_shape.first * _shape.second];
		return arr[i * _shape.second + j];
	}
	const T &operator()(size_t i, size_t j) const {
		if (i > _shape.first || j > _shape.second)
			return arr[_shape.first * _shape.second];
		return arr[i * _shape.second + j];
	}
};

#endif // FOCUS_INCLUDES_TYPES_MATRIX_2D_HPP

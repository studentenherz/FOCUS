#if !defined(FOCUS_INCLUDE_TYPES_VECTOR_HPP)
#define FOCUS_INCLUDE_TYPES_VECTOR_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>

/**
 * n-Vector class
 * Fixed size array
 */
template<size_t n>
class Vector{
	double v[n + 1];
public:
	/**
	 * Default constructor creates the zero vector
	 */
	Vector() {
		for (size_t i = 0; i < n; i++)
			v[i] = 0;
	}

	// /**
	//  * Create a vector from x array
	//  */
	// Vector(double x, ...){
	// 	std::va_list args;
	// 	va_start(args, x);
	// 	v[0] = x;
	// 	for (size_t i = 1; i < n; i++)
	// 		v[i] = va_arg(args, double);
	// 	va_end(args);
	// }

	Vector(std::initializer_list<double> l){
		size_t index = 0;
		for(const double* it = l.begin(); it != l.end(); it++)
			v[index++] = *it;
	}

	/**
	 * Access to vector component
	 */
	double &operator[](size_t i) {return v[(i < n ? i : n)];}
	const double &operator[](size_t i) const {return v[(i < n ? i : n)];}
};

// Vector Algebra

/**
 * Vectorial sum
 * @param a one vector
 * @param b other vector
 * @return vectorial sum a + b
 */
template<size_t n>
Vector<n> operator+(Vector<n> a, Vector<n> b){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] + b[i];
	return c;
}

/**
 * Vectorial difference
 * @param a one vector
 * @param b another vector
 * @return vectorial sum a + (-b)
 */
template<size_t n>
Vector<n> operator-(Vector<n> a, Vector<n> b){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] - b[i];
	return c;
}

/**
 * Scalar vector multiplication
 * @param a a vector
 * @param t a scalar
 * @return scaled vector t * a
 */
template<size_t n>
Vector<n> operator*(Vector<n> a, double t){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] * t;
	return c;
}

/**
 * Scalar vector multiplication
 * @param t a scalar
 * @param a a vector
 * @return scaled vector t * a
 */
template<size_t n>
Vector<n> operator*(double t, Vector<n> a){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] * t;
	return c;
}

/**
 * Vector divided by scalar
 * @param a a vector
 * @param t a scalar
 * @return scaled vector a / t
 */
template<size_t n>
Vector<n> operator/(Vector<n> a, double t){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] / t;
	return c;
}

/**
 * Dot product between two vector
 * @param a one vector
 * @param b another vector
 * @return dot product a . b
 */
template<size_t n>
double dot(Vector<n> a, Vector<n> b){
	double s = 0;
	for (size_t i = 0; i < n; i++)
		s += a[i] * b[i];
	return s;
}


/**
 * Modulus of vector (2-norm)
 * @param v a vector
 * @return |v|
 */
template<size_t n>
double mod(Vector<n> v){
	return sqrt(dot(v, v));
}

// IO for vectors

/**
 * Out stream operator
 */
template<size_t n>
std::ostream& operator<<(std::ostream& out, const Vector<n>& v){
	for(size_t i = 0; i < n; i++)
		out << v[i] << ' ';
	return out; 
}

/**
 * In stream operator
 */
template<size_t n>
std::istream& operator>>(std::istream& in, Vector<n>& v){
	for(size_t i = 0; i < n; i++)
		in >> v[i];
	return in; 
}


// 3-Vector
typedef Vector<3> Vector3;

/**
 * Cross product between two 3-Vectors
 * @param a one vector
 * @param b another vector
 * @return dot product a x b
 */
Vector3 cross(Vector3 a, Vector3 b){
	Vector3 c;
	for (size_t i = 0; i < 3; i++)
		c[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
	return c;
}

// State: space + velocity 3-Vectors
typedef Vector<6> State;

/**
 * Get position form a state
 * @param x state
 * @return position
 */
Vector3 get_position(State x){
	Vector3 r;
	for(size_t i = 0; i < 3; i++)
		r[i] = x[i];
	return r;
}

/**
 * Get velocity form a state
 * @param x state
 * @return velocity
 */
Vector3 get_velocity(State x){
	Vector3 r;
	for(size_t i = 3; i < 6; i++)
		r[i - 3] = x[i];
	return r;
}

#endif // FOCUS_INCLUDE_TYPES_VECTOR_HPP

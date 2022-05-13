#if !defined(FOCUS_INCLUDE_TYPES_VECTOR_HPP)
#define FOCUS_INCLUDE_TYPES_VECTOR_HPP

#include <iostream>
#include <cmath>

/**
 * 3-Vector class with standard mathematical operaions
 */
class Vector{
	double v[4];
public:
	/**
	 * Default constructor creates the zero vector
	 */
	Vector(): v{0, 0, 0} {}

	/**
	 * Create a vector (x, y, z)
	 */
	Vector(double x, double y, double z): v{x, y, z} {}

	/**
	 * Access to vector component
	 */
	double &operator[](size_t i) {return v[(i < 3 ? i : 3)];}
	const double &operator[](size_t i) const {return v[(i < 3 ? i : 3)];}
};

// Vector Algebra

/**
 * Vectorial sum
 * @param a one vector
 * @param b other vector
 * @return vectorial sum a + b
 */
Vector operator+(Vector a, Vector b){
	Vector c;
	for(size_t i = 0; i < 3; i++)
		c[i] = a[i] + b[i];
	return c;
}

/**
 * Vectorial difference
 * @param a one vector
 * @param b another vector
 * @return vectorial sum a + (-b)
 */
Vector operator-(Vector a, Vector b){
	Vector c;
	for(size_t i = 0; i < 3; i++)
		c[i] = a[i] - b[i];
	return c;
}

/**
 * Scalar vector multiplication
 * @param a a vector
 * @param t a scalar
 * @return scaled vector t * a
 */
Vector operator*(Vector a, double t){
	Vector c;
	for(size_t i = 0; i < 3; i++)
		c[i] = a[i] * t;
	return c;
}

/**
 * Scalar vector multiplication
 * @param t a scalar
 * @param a a vector
 * @return scaled vector t * a
 */
Vector operator*(double t, Vector a){
	Vector c;
	for(size_t i = 0; i < 3; i++)
		c[i] = a[i] * t;
	return c;
}

/**
 * Vector divided by scalar
 * @param a a vector
 * @param t a scalar
 * @return scaled vector a / t
 */
Vector operator/(Vector a, double t){
	Vector c;
	for(size_t i = 0; i < 3; i++)
		c[i] = a[i] / t;
	return c;
}

/**
 * Dot product between two vector
 * @param a one vector
 * @param b another vector
 * @return dot product a . b
 */
double dot(Vector a, Vector b){
	double s = 0;
	for (size_t i = 0; i < 3; i++)
		s += a[i] * b[i];
	return s;
}

/**
 * Cross product between two vector
 * @param a one vector
 * @param b another vector
 * @return dot product a x b
 */
Vector cross(Vector a, Vector b){
	Vector c;
	for (size_t i = 0; i < 3; i++)
		c[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
	return c;
}

/**
 * Modulus of vector (2-norm)
 * @param v a vector
 * @return |v|
 */
double mod(Vector v){
	return sqrt(dot(v, v));
}

// IO for vectors

/**
 * Out stream operator
 */
std::ostream& operator<<(std::ostream& out, const Vector& v){
	for(size_t i = 0; i < 3; i++)
		out << v[i] << ' ';
	return out; 
}

/**
 * In stream operator
 */
std::istream& operator>>(std::istream& in, Vector& v){
	for(size_t i = 0; i < 3; i++)
		in >> v[i];
	return in; 
}


#endif // FOCUS_INCLUDE_TYPES_VECTOR_HPP

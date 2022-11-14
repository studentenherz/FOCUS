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
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Vector() {
		for (size_t i = 0; i < n; i++)
			v[i] = 0;
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Vector(std::initializer_list<double> l){
		size_t index = 0;
		for(const double* it = l.begin(); it != l.end(); it++)
			v[index++] = *it;
	}

	/**
	 * Access to vector component
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double &operator[](size_t i) {return v[(i < n ? i : n)];}
	
	#ifdef __CUDACC__
	__host__ __device__
	#endif
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator+(const Vector<n>& a, const Vector<n>& b){
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator-(const Vector<n>& a, const Vector<n>& b){
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator*(const Vector<n>& a, double t){
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator*(double t, const Vector<n>& a){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] * t;
	return c;
}

/**
 * Element-wise vector multiplication
 * @param a a vector
 * @param b other vector
 * @return scaled vector t * a
 */
template<size_t n>
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator*(const Vector<n>& a, const Vector<n>& b){
	Vector<n> c;
	for(size_t i = 0; i < n; i++)
		c[i] = a[i] * b[i];
	return c;
}

/**
 * Vector divided by scalar
 * @param a a vector
 * @param t a scalar
 * @return scaled vector a / t
 */
template<size_t n>
#ifdef __CUDACC__
__host__ __device__
#endif
Vector<n> operator/(const Vector<n>& a, double t){
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
#ifdef __CUDACC__
__host__ __device__
#endif
double dot(const Vector<n>& a, const Vector<n>& b){
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
#ifdef __CUDACC__
__host__ __device__
#endif
double mod(const Vector<n>& v){
	return sqrt(dot(v, v));
}

/**
 * Checks if there is a nan value inside
 * @param v Vector
 * @return true if there is a nan value
 */
template<size_t n>
#ifdef __CUDACC__
__host__ __device__
#endif
bool hasnan(const Vector<n>& v){
	for(size_t i = 0; i<n; i++)
		if(std::isnan(v[i]))
			return true;
	return false;
}

// IO for vectors

/**
 * Out stream operator
 */
template<size_t n>
#ifdef __CUDACC__
__host__
#endif 
std::ostream& operator<<(std::ostream& out, const Vector<n>& v){
	for(size_t i = 0; i < n; i++)
		out << v[i] << ' ';
	return out; 
}

/**
 * In stream operator
 */
template<size_t n>
#ifdef __CUDACC__
__host__
#endif 
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector3 cross(const Vector3& a, const Vector3& b){
	Vector3 c;
	for (size_t i = 0; i < 3; i++)
		c[i] = a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
	return c;
}

/**
 * Angle between two 3-Vectors
 * @param a one vector
 * @param b another vector
 * @return dot product a x b
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double angle_between(const Vector3& a, const Vector3& b){
	return acos(dot(a, b) / (mod(a) * mod(b)));
}

/**
 * Pitch between two 3-Vectors
 * @param a one vector
 * @param b another vector
 * @return dot product a x b
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double pitch_between(const Vector3& a, const Vector3& b){
	return dot(a, b) / (mod(a) * mod(b));
}

// State: space + velocity 3-Vectors
typedef Vector<6> State;

/**
 * Get position form a state
 * @param x state
 * @return position
 */
#ifdef __CUDACC__
__host__ __device__
#endif
Vector3 get_position(const State& x){
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
#ifdef __CUDACC__
__host__ __device__
#endif
Vector3 get_velocity(const State& x){
	Vector3 r;
	for(size_t i = 3; i < 6; i++)
		r[i - 3] = x[i];
	return r;
}

/**
 * Convert vector from cylindrical coordinates to cartesian
 * @param v cylindrical vector
 * @param r cylindrical position
 * @return vector in cartesian
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline Vector3 cyl2cart(const Vector3& v, double theta){
	Vector3 cart;
	double c = cos(theta), s = sin(theta);

	cart[0] = v[0] * c - v[1] * s;
	cart[1] = v[0] * s + v[1] * c;
	cart[2] = v[2];

	return cart;
}

/**
 * Convert vector from cartesian coordinates to cylindrical
 * @param v cartesian vector
 * @param r cylindrical position
 * @return vector in cylindrical
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline Vector3 cart2cyl(const Vector3& v, double theta){
	Vector3 cyl;
	double c = cos(theta), s = sin(theta);

	cyl[0] = v[0] * c + v[1] * s;
	cyl[1] = v[1] * c - v[0] * s;
	cyl[2] = v[2];

	return cyl;
}

/**
 * Convert from cylindrical coordinates to cartesian
 * @param r cylindrical vcoordinates
 * @return cartesian coordinates
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline Vector3 cyl2cart(const Vector3& r){
	Vector3 cart;

	cart[0] = r[0] * cos(r[1]);
	cart[1] = r[0] * sin(r[1]);
	cart[2] = r[2];

	return cart;
}

/**
 * Convert from cartesian coordinates to cylindrical
 * @param r cartesian coordinates
 * @return cylindrical coordinates
 */
#ifdef __CUDACC__
__host__ __device__
#endif
inline Vector3 cart2cyl(const Vector3& r){
	Vector3 cyl;

	cyl[0] = sqrt(r[0] * r[0] + r[1] * r[1]);
	cyl[1] = atan2(r[1], r[0]);
	cyl[2] = r[2];

	return cyl;
}


#endif // FOCUS_INCLUDE_TYPES_VECTOR_HPP

#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "interpolations.hpp"
#include "types/vector.hpp"
#include "types/array.hpp"

struct Particle{
	double m;			///< mass
	int q;				///< charge
	ulong n;			///< Principal quantum number
	double t;			///< Time since it's been in this quantum state

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Particle() {}

	/**
	 * Constructor
	 * @param q charge
	 * @param m mass
	 * @param n Principal quantum number
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Particle(double m, uint q, uint n = 1): m(m), q(q), n(n) {}
};

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

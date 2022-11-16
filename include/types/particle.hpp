#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "interpolations.hpp"
#include "types/vector.hpp"
#include "types/array.hpp"

struct Particle{
	uint q;				///< charge
	double m;			///< mass
	uint n;				///< Principal quantum number

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
	Particle(uint q, double m, uint n = 0): q(q), m(m), n(n) {}
};

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

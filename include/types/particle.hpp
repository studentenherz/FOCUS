#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "interpolations.hpp"
#include "types/vector.hpp"
#include "types/array.hpp"

struct Particle{
	double q;			///< charge
	double m;			///< mass

	Particle() {}

	/**
	 * Constructor
	 * @param q charge
	 * @param m mass
	 */
	Particle(double q, double m): q(q), m(m) {}
};

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

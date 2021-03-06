#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "types/vector.hpp"

struct ParticleSpecies{
	double q;			///< charge
	double m;			///< mass
	double logl;	///< Lorentz logarithm

	double (*T)(Vector3, double); ///< temperature profile
	double (*n)(Vector3, double); ///< density profile

	ParticleSpecies() {}

	ParticleSpecies(double charge, double mass, double lorentz_logarithm, double (*Temperature)(Vector3, double), double (*density)(Vector3, double)): q(charge), m(mass), logl(lorentz_logarithm), T(Temperature), n(density) {}
};


#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

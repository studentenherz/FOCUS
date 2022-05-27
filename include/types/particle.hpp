#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "types/vector.hpp"

struct ParticleSpecies{
	double q;			// charge
	double m;			// mass
	double logl;	// Lorentz logarithm

	double (*T)(Vector3, double); // temperature profile
	double (*n)(Vector3, double); // density profile
};


#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

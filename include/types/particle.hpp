#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "types/vector.hpp"

/**
 * Base particle species class
 */
struct ParticleSpecies{
	double q;			///< charge
	double m;			///< mass
	double logl;	///< Lorentz logarithm

	ParticleSpecies() {}

	/**
	 * Constructor
	 * @param charge charge of the particles
	 * @param mass mass of the particles
	 * @param lorentz_logarith Lorentz's logarithm of the species
	 */
	ParticleSpecies(double charge, double mass, double lorentz_logarithm): q(charge), m(mass), logl(lorentz_logarithm) {}
	virtual ~ParticleSpecies() = default;
	
	virtual double T(Vector3, double) = 0; ///< temperature profile
	virtual double n(Vector3, double) = 0; ///< density profile

};


/**
 * Particle species class that implements
 * a constant temperature and density profile.
 */
class ConstProfileParticle: public ParticleSpecies{
	double T_;
	double n_;
public:
	/**
	 * Constructor
	 * @param charge charge of the particles
	 * @param mass mass of the particles
	 * @param lorentz_logarith Lorentz's logarithm of the species
	 * @param temperature constant temperature of the particle ensemble
	 * @param density constant density profile of the species
	 */
	ConstProfileParticle(double charge, double mass, double lorentz_logarithm, double temperature, double density) : ParticleSpecies(charge, mass, lorentz_logarithm), T_(temperature), n_(density) {}

	double T(Vector3 /* r */, double /* t */){
		return T_;
	}

	double n(Vector3 /* r */, double /* t */){
		return n_;
	}
};

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "interpolations.hpp"
#include "types/vector.hpp"
#include "types/array.hpp"

/**
 * Base particle species class
 */
struct ParticleSpecies{
	double q;			///< charge
	double m;			///< mass
	double logl;	///< Lorentz logarithm

	Array<double> psi_prof;
	Array<double> T_prof;
	Array<double> n_prof;

	ParticleSpecies() {}

	/**
	 * Constructor
	 * @param charge charge of the particles
	 * @param mass mass of the particles
	 * @param lorentz_logarith Lorentz's logarithm of the species
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	ParticleSpecies(double charge, double mass, double lorentz_logarithm, Array<double> psi, Array<double> Temp, Array<double> density): q(charge), m(mass), logl(lorentz_logarithm), psi_prof(psi), T_prof(Temp), n_prof(density) {}

	/**
	 * Construct in host for device 
	 * @param other 
	 */
	#ifdef __CUDACC__
	__host__
	void construct_in_host_for_device( Array<double> psi, Array<double> Temp, Array<double> density){
		psi_prof.construct_in_host_for_device(psi);
		T_prof.construct_in_host_for_device(Temp);
		n_prof.construct_in_host_for_device(density);
	}
	#endif
	
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double T(double psi, double /* t */){
		if (psi_prof.size() < 3)
			return T_prof[0];
		return lagrange_interpolation_3(psi, psi_prof, T_prof);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double n(double psi, double /* t */){
		if (psi_prof.size() < 3)
			return n_prof[0];
		return lagrange_interpolation_3(psi, psi_prof, n_prof);
	}

};


// /**
//  * Particle species class that implements
//  * a constant temperature and density profile.
//  */
// class ConstProfileParticle: public ParticleSpecies{
// 	double T_;
// 	double n_;
// public:
// 	/**
// 	 * Constructor
// 	 * @param charge charge of the particles
// 	 * @param mass mass of the particles
// 	 * @param lorentz_logarith Lorentz's logarithm of the species
// 	 * @param temperature constant temperature of the particle ensemble
// 	 * @param density constant density profile of the species
// 	 */
// 	#ifdef __CUDACC__
// 	__host__ __device__
// 	#endif
// 	ConstProfileParticle(double charge, double mass, double lorentz_logarithm, double temperature, double density) : ParticleSpecies(charge, mass, lorentz_logarithm), T_(temperature), n_(density) {}

// 	#ifdef __CUDACC__
// 	__host__ __device__
// 	#endif
// 	double T(Vector3 /* r */, double /* t */){
// 		return T_;
// 	}

// 	#ifdef __CUDACC__
// 	__host__ __device__
// 	#endif
// 	double n(Vector3 /* r */, double /* t */){
// 		return n_;
// 	}
// };

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

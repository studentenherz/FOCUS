#if !defined(FOCUS_INCLUDE_TYPES_PARTICLE_HPP)
#define FOCUS_INCLUDE_TYPES_PARTICLE_HPP

#include "interpolations.hpp"
#include "types/vector.hpp"
#include "types/array.hpp"

double logl_prefactor = 18.4527;

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

/**
 * Base particle species class
 */
struct PlasmaParticleSpecies: public Particle{

	Array<double>* psi_prof;
	Array<double>* T_prof;
	Array<double>* n_prof;

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	PlasmaParticleSpecies() {}

	/**
	 * Constructor
	 * @param q charge of the particles
	 * @param m mass of the particles
	 * @param psi Poloidal flux profile
	 * @param Temp Temperature profile [keV]
	 * @param density Particle density [10^19/m^3]
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	PlasmaParticleSpecies(double q, double m, Array<double>& psi, Array<double>& Temp, Array<double>& density): Particle(q, m), psi_prof(&psi), T_prof(&Temp), n_prof(&density) {}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double T(double psi, double /* t */){
		if (psi_prof->size() < 3)
			return T_prof->operator[](0);
		return lagrange_interpolation_3(psi, *psi_prof, *T_prof);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double n(double psi, double /* t */){
		if (psi_prof->size() < 3)
			return n_prof->operator[](0);
		return lagrange_interpolation_3(psi, *psi_prof, *n_prof);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double logl(Particle a, double psi, double t){
		return logl_prefactor + log(T(psi, t) / (n(psi, t) * sqr(1 + m/a.m) * sqr(a.q) * sqr(sqr(q))));
	}

};

#endif // FOCUS_INCLUDE_TYPES_PARTICLE_HPP

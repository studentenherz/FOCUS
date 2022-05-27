#if !defined(FOCUS_INCLUDE_COLLISIONS_HPP)
#define FOCUS_INCLUDE_COLLISIONS_HPP

#include <cmath>

#include "util.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "random.hpp"

/**
 * Error function minus it's derivative
 */
double erf_minus_d_erf(double x){
	return (erf(x) - x * two_over_sqrt_pi * exp(-sqr(x)));
}

/**
 * Another definition for formulae simplicity
 */
double G(double x){
	return erf_minus_d_erf(x)/(2 * sqr(x));
}

/**
 * This class implements the elastic collisions
 * calculation from Focker-Planks' theory in the
 * Langevin equation using Îto's method.
 */
class FockerPlank{
	Array<ParticleSpecies> beta;	// Particle species involved
	ParticleSpecies alpha;				// Test particle species

	double _eta;				// Dimensionless constant
	NormalRand gauss;	// Gaussian random generator
public:
	FockerPlank(unsigned long long seed, Array<ParticleSpecies> plasma_particles, ParticleSpecies test_particle, double eta): gauss(seed), beta(plasma_particles), alpha(test_particle), _eta(eta) {}

	/**
	 * Slowing down from elastic collisions
	 * @param x current state
	 * @param t current time
	 * @return rate of change of the parallel velocity
	 */
	Vector3 slow_down(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		double nu_sd = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i].T(r, t);
			nu_sd +=  sqr(beta[i].q) *  beta[i].n(r, t) * (1 + alpha.m/beta[i].m) * beta[i].logl * erf_minus_d_erf(xb);
		}

		// other terms
		nu_sd *= _eta * sqr(alpha.q) / (pow(v_mod, 3) * sqr(alpha.m));

		return - nu_sd * v;
	}

	/**
	 * Coefficient of dispersion in the parallel direction
	 * @param x current state
	 * @param t current time
	 * @return parallel dispersion coefficient
	 */
	double parallel_dispersion_coeff(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i].T(r, t);
			nu +=  sqr(beta[i].q) * beta[i].n(r, t) * beta[i].logl * G(xb);
		}

		// other terms
		nu *= 2 * _eta * sqr(alpha.q) / (pow(v_mod, 3) * sqr(alpha.m));
		
		return nu;
	}

	/**
	 * Coefficient of dispersion in the perpendicular direction
	 * @param x current state
	 * @param t current time
	 * @return perpendicular dispersion coefficient
	 */
	double perpendicular_dispersion_coeff(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i].T(r, t);
			nu +=  sqr(beta[i].q) * beta[i].n(r, t) * beta[i].logl * (erf(xb) - G(xb));
		}

		// other terms
		nu *= 2 * _eta * sqr(alpha.q) / (pow(v_mod, 3) * sqr(alpha.m));
		
		return nu;
	}

	/**
	 * Stochastic dispersion from elastic collisions
	 * @param x current state
	 * @param t current time
	 * @return rate of dispersion of the velocity
	 */ 
	Vector3 dispersion(const State& x, double t){
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		// Dispersion coefficients
		double nu_ll = parallel_dispersion_coeff(x, t);
		double nu_T = perpendicular_dispersion_coeff(x, t);

		// Parallel (ll) unit vector
		Vector3 e_ll = v / v_mod;

		// First perpendicular (T) unit vector
		Vector3 e_x(1.0, 0.0, 0.0);
		Vector3 e_T_1 = cross(e_ll, e_x);
		if (mod(e_T_1) == 0){
			Vector3 e_y(0.0, 1.0, 0.0);
			e_T_1 = cross(e_ll, e_y);
		}
		e_T_1 = e_T_1 / mod(e_T_1);

		// Second perpendicular (T) unit vector
		Vector3 e_T_2 = cross(e_ll, e_T_1);

		Vector3 dvdt = sqrt(nu_ll) * v_mod * gauss() * e_ll + sqrt(nu_T / 2) * v_mod * (gauss() * e_T_1 + gauss() * e_T_2);

		return dvdt;
	}

	/**
	 * Euler stepper from the Îto calculation
	 * @param[out] x current state
	 * @param t current time
	 * @param dt time delta
	 */
	void euler_step(State& x, double t, double dt){
		Vector3 dv = slow_down(x, t) * dt + dispersion(x, t) * sqrt(dt);
	
		x[3] += dv[0];
		x[4] += dv[1];
		x[5] += dv[2];
	}
};

#endif // FOCUS_INCLUDE_COLLISIONS_HPP

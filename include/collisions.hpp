/**
 * @file collisions.hpp
 * @brief Implementation of the Focker-Plank elastic collisions theory.
 * 
 * The Focker-Plank theory [1] describes the slowing down and dispersion of ions in a plasma
 * taking in consideration that, given the conditions of a fusion plasma, it is far more
 * probable that a large deviation is cause by multiple succesive small deviation rather
 * than a large unique one.
 * 
 * Using the Itô's calculus one gets that the variation of the velocity \f$ v \f$ in the 
 * direction \f$ i \f$ caused by the elastic collisions is
 * 
 * \f[
 * 	\frac{dv_i}{dt} = F_i(v, t) + \sqrt{D_{ii}(v, t)} \xi_i(t)
 * \f]
 * 
 * where \f$ \xi_i(t) \f$ is white noise with Gaussian distribution that fulfills \f$ \langle \xi_i(t) \rangle = 0 \f$ 
 * and \f$ \langle \xi_i(t) \xi_k(t') \rangle = \delta(t-t')\delta_{ik} \f$. For a plasma with species \f$ \beta \f$
 * in a maxwellian equilibrium the friction coefficient \f$ F_i \f$ and the difusion tensor \f$ D_{ii} \f$ are known:
 * 
 * \f[
 * 	\begin{aligned}
 *	 F_{||}(v) &= -\nu_\text{sd}(v) v\\
 * 	 D_{||}(v) &= \nu_{||}(v) v^2 \\
 *   D_\perp(v) &= \nu_\perp(v) v^2
 * \end{aligned}
 * \f]
 * 
 * where the slowind down factor is given by
 * 
 * \f[
 * 	\nu_\text{sd}(v) = \sum\limits_\beta \frac{A_\text{D}^\beta}{ 2 v^3} \left(1+ 
 * 		\frac{m_\alpha}{m_\beta}\right)\left(\phi(x_\beta) - x_\beta \phi'(x_\beta)\right)
 * \f]
 * 
 * and the dispersion frequencies along the instant parallel and perpendicular directions 
 * are given respectively by
 * 
 * \f[
 * 	  \nu_{||}(v) = \sum\limits_\beta \frac{A_\text{D}^\beta}{v^3}G(x_\beta)
 * \f]
 * 
 * and
 * 
 * \f[
 * 	  \nu_\perp(v) = \sum\limits_\beta \frac{A_\text{D}^\beta}{v^3}\left(\phi(x_\beta) - G(x_\beta)\right).
 * \f]
 * 
 * In the previous equations \f$ x_\beta = v/v_{s,\beta} \f$ where \f$ v_{s,\beta} = \sqrt{2 k_\text{B} T / m_\beta} \f$,
 *  \f$ \phi(x) \f$ is the error function, 
 * 
 * \f[
 * 	G(x) = \frac{\phi(x) - x \phi'(x)}{2x^2}
 * \f]
 * 
 * and 
 * 
 * \f[
 *	A_\text{D}^\beta = \frac{q_\alpha^2q_\beta^2}{2 \pi \epsilon_0^2 m _\alpha^2} n_\beta \ln \Lambda_\beta.
 * \f]
 * 
 * 
 * [1] Krall, N., Trivelpiece, A., Kempton, J. Principles of Plasma Physics. International
 * series in pure and applied physics. McGraw-Hill, 1973. URL https://books.google.com.ar/books?id=b0BRAAAAMAAJ.
 */ 

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
	Array<ParticleSpecies*> beta;	// Particle species involved
	ParticleSpecies& alpha;				// Test particle species

	double _eta;				// Dimensionless constant
	NormalRand gauss;	// Gaussian random generator
public:
	FockerPlank(unsigned long long seed, Array<ParticleSpecies*> plasma_particles, ParticleSpecies& test_particle, double eta): beta(plasma_particles), alpha(test_particle), _eta(eta), gauss(seed) {}

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
			double xb = v_mod / beta[i]->T(r, t);
			nu_sd +=  sqr(beta[i]->q) *  beta[i]->n(r, t) * (1 + alpha.m/beta[i]->m) * beta[i]->logl * erf_minus_d_erf(xb);
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
			double xb = v_mod / beta[i]->T(r, t);
			nu +=  sqr(beta[i]->q) * beta[i]->n(r, t) * beta[i]->logl * G(xb);
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
			double xb = v_mod / beta[i]->T(r, t);
			nu +=  sqr(beta[i]->q) * beta[i]->n(r, t) * beta[i]->logl * (erf(xb) - G(xb));
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
		Vector3 e_x = {1, 0, 0};
		Vector3 e_T_1 = cross(e_ll, e_x);
		if (mod(e_T_1) == 0){
			Vector3 e_y = {0, 1, 0};
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

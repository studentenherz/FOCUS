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
 * With this we can include the effect of collisions for a single particle using Langevin equations
 * 
 * \f[
 * 	\begin{aligned}
 *		\frac{dv_1}{dt} &= -\nu_\text{sd} v + \sqrt{\nu_{||} v^2} \xi_1(t) \\ 
 * 		\frac{dv_{2, 3}}{dt} &= \sqrt{\frac{\nu_\perp}{2} v^2} \xi_{2, 3}(t),
 * 	\end{aligned}
 * \f]
 * 
 * where 1 represents the direction instantly parallel to the velocity and 2, 3 represent a pair of perpendicular
 * directions that make a right handed coordinates system. In the doctoral thesis by Clauser [2] it is shown
 * that these equations can be solved using an Euler step wihout compromising the accuracy of the calculation. Hence
 * the variations are calculated using
 * 
 * \f[
 * 	\begin{aligned}
 * 		\Delta v_1 &= -\nu_\text{sd} v \Delta t_\text{col} + \sqrt{\nu_{||} \Delta t_\text{col}} v N_1\\
 * 		\Delta v_{2, 3} &= \sqrt{\frac{\nu_\perp \Delta t_\text{col}}{2} } v N_{2, 3}
 * 	\end{aligned}
 * \f]
 * 
 * where \f$ N_i \f$ are random numbers with a Gaussian distribution with mean 0 and variance 1.
 * 
 * As in other parts of the code, for calculations purposes it is wise to use some dimensionless factor to prevent
 * floating point error, in this case we have to calculate \f$ \nu \Delta t \f$ that is a dimensionless quantity
 * 
 * \f[
 * 	\begin{aligned}
 * 		\nu \Delta t &= \frac{A^\beta_D}{2 v^3} F(x_\beta) \Delta t \\
 * 		             &= \frac{q_\alpha^2q_\beta^2 \Delta t}{4 \pi \epsilon_0^2 m _\alpha^2 v^3} n_\beta \ln \Lambda_\beta F(x_\beta) \\
 * 		             &= \frac{e^4Z_\alpha^2Z_\beta^2 \tau \Delta \hat{t}}{4 \pi \epsilon_0^2 m_e^2 \hat{m} _\alpha^2 v_0^3 \hat{v}^3} n_0 \hat{n}_\beta \ln \Lambda_\beta F(x_\beta) \\
 *		             &= \eta \frac{Z_\alpha^2Z_\beta^2\Delta\hat{t}}{\hat{m} _\alpha^2\hat{v}^3} \hat{n}_\beta \ln \Lambda_\beta F(x_\beta)
 * 	\end{aligned}
 * \f]
 * 
 * where \f$ F(x_\beta) \f$ depends on the component we are calculating and defining the dimensionless factor
 * 
 * \f[
 * 	\eta = \frac{e^4 \tau n_0}{4 \pi \epsilon_0^2 m_e^2 v_0^3}.
 * \f]
 * 
 * [1] Krall, N., Trivelpiece, A., Kempton, J. Principles of Plasma Physics. International
 * series in pure and applied physics. McGraw-Hill, 1973. URL https://books.google.com.ar/books?id=b0BRAAAAMAAJ.
 * 
 * [2] Clauser, C. F. Dinámica de partículas alfa en plasmas magnetizados y el efecto
 * de las colisiones en la interacción partícula-plasma. Tesis Doctoral, 2018.
 */ 

#if !defined(FOCUS_INCLUDE_COLLISIONS_HPP)
#define FOCUS_INCLUDE_COLLISIONS_HPP

#include <cmath>

#include "util.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "random.hpp"
#include "magnetic_field.hpp"

/**
 * Error function minus it's derivative
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double erf_minus_d_erf(double x){
	return (erf(x) - x * two_over_sqrt_pi * exp(-sqr(x)));
}

/**
 * Another definition for formulae simplicity
 */
#ifdef __CUDACC__
__host__ __device__
#endif
double G(double x){
	return erf_minus_d_erf(x)/(2 * sqr(x));
}

/**
 * This class implements the elastic collisions
 * calculation from Focker-Planks' theory in the
 * Langevin equation using Îto's method.
 */
template<typename NormalRand_t>
class FockerPlank{
	Array<ParticleSpecies*>& beta;	// Particle species involved
	ParticleSpecies& alpha;				// Test particle species

	MagneticFieldFromMatrix& _B;

	double _eta;				// Dimensionless constant
	NormalRand_t& gauss;	// Gaussian random generator
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	FockerPlank(Array<ParticleSpecies*>& plasma_particles, ParticleSpecies& test_particle, MagneticFieldFromMatrix& B, double eta, NormalRand_t& normal_rand): beta(plasma_particles), alpha(test_particle), _B(B), _eta(eta), gauss(normal_rand) {}

	/**
	 * Slowing down from elastic collisions
	 * @param x current state
	 * @param t current time
	 * @return rate of change of the parallel velocity
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Vector3 slow_down(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		double psi = _B.psi(r, t);

		double nu_sd = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i]->T(psi, t);
			nu_sd +=  sqr(beta[i]->q) *  beta[i]->n(psi, t) * (1 + alpha.m/beta[i]->m) * beta[i]->logl * erf_minus_d_erf(xb);
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
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double parallel_dispersion_coeff(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);

		double psi = _B.psi(r, t);

		double nu = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i]->T(psi, t);
			nu +=  sqr(beta[i]->q) * beta[i]->n(psi, t) * beta[i]->logl * G(xb);
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
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double perpendicular_dispersion_coeff(const State& x, double t){
		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);
		double v_mod = mod(v);
		
		double psi = _B.psi(r, t);

		double nu = 0;
		// terms that depend on plasma particles
		for(size_t i = 0; i < beta.size(); i++){
			double xb = v_mod / beta[i]->T(psi, t);
			nu +=  sqr(beta[i]->q) * beta[i]->n(psi, t) * beta[i]->logl * (erf(xb) - G(xb));
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
	#ifdef __CUDACC__
	__host__ __device__
	#endif
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
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void euler_step(State& x, double t, double dt){
		Vector3 dv = slow_down(x, t) * dt + dispersion(x, t) * sqrt(dt);
	
		x[3] += dv[0];
		x[4] += dv[1];
		x[5] += dv[2];
	}
};

#endif // FOCUS_INCLUDE_COLLISIONS_HPP

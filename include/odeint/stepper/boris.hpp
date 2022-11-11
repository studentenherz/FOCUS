#if !defined(FOCUS_INCLUDE_ODEINT_STEPPER_BORIS_HPP)
#define FOCUS_INCLUDE_ODEINT_STEPPER_BORIS_HPP

#include "types/vector.hpp"
#include "lorentz.hpp"

/**
 * Boris algorithm for integrating charged particles subjected to Lorentz force [1, 2]
 * 
 * [1] [https://www.particleincell.com/2011/vxb-rotation/](https://www.particleincell.com/2011/vxb-rotation/)
 * 
 * [2] Qin, H., Zhang, S., Xiao, J., Liu, J., Sun, Y., & Tang, W. M. (2013). Why is Boris algorithm so good?. Physics of Plasmas, 20(8), 084503. [https://doi.org/10.1063/1.4818428](https://doi.org/10.1063/1.4818428)
 */
template<typename system_type, typename state_type = State, typename scalar_type = double>
class Boris{
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Boris() {}
	
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void do_step(system_type sys, state_type& x, scalar_type time, scalar_type dt){
		Vector3 r_cyl = get_position(x);
		Vector3 v_cyl = get_velocity(x);

		// Convert to cartesian
		Vector3 rc = cyl2cart(r_cyl);
		Vector3 vc = cyl2cart(v_cyl, r_cyl[1]);

		/* 
			The Boris method is similar to a leapfrog method as is has
			staggered position and velocities. Some define then the velocity
			v_k as v(t_k - dt/2) and position x_k as x(t_k). In order to be coherent
			with our code and have v_k to be v(t_k) we just need to be careful and
			and calculate the fields at t_k + dt/2.
		*/

		// Update position by dt / 2
		rc = rc + sys.gam * dt * vc / 2;
		r_cyl = cart2cyl(rc);

		// Calculate the fields at x(t + dt)
		Vector3 Bcil = sys.B(r_cyl, time + dt / 2);
		Vector3 Ecil = sys.E(r_cyl, time + dt / 2);

		// Convert to cartesian
		Vector3 Bc = cyl2cart(Bcil, r_cyl[1]);
		Vector3 Ec = cyl2cart(Ecil, r_cyl[1]);

		// Calculate v'
		Vector3 vp = vc + sys.gam * dt * (Ec  + cross(vc, Bc) / 2);
		Vector3 t = sys.Z_m * Bc * dt / 2;

		double s = 1.0/(1.0 + dot(t, t));
		double v_dot_t = vp[0]*t[0] + vp[1]*t[1] + vp[2]*t[2];

		// Calculate v(t + dt)
		vc = s * (vp + v_dot_t * t + cross(vp, t));

		// Update position by dt/2 with new velocity
		rc = rc + sys.gam * vc * dt / 2;

		// Convert back to cyl
		r_cyl = cart2cyl(rc);
		v_cyl = cart2cyl(vc, r_cyl[1]);

		// Update particle state
		x[0] = r_cyl[0];
		x[1] = r_cyl[1];
		x[2] = r_cyl[2];
		x[3] = v_cyl[0];
		x[4] = v_cyl[1];
		x[5] = v_cyl[2];
	}
};

#endif // FOCUS_INCLUDE_ODEINT_STEPPER_BORIS_HPP

#if !defined(FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP)
#define FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP

/**
 * Euler integrator order 1. Just for testing purposes
 * as it is quite unprecise.
 */
template<typename system_type, typename state_type, typename scalar_type>
class EulerStepper{
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	EulerStepper() {}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	static void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		state_type dxdt;
		sys(x, dxdt, t); // dxdt = f(x, t)

		x = x + dxdt * dt;
	}
};

#endif // FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP

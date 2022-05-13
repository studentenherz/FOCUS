#if !defined(FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP)
#define FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP

/**
 * Euler integrator order 1. Just for testing purposes
 * as it is quite unprecise.
 */
template<typename system_type, typename state_type, typename derivative_type, typename scalar_type>
class EulerStepper{
public:
	EulerStepper() {}
	static void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		derivative_type dxdt;
		sys(x, dxdt, t); // dxdt = f(x, t)

		x = x + dxdt * dt;
	}
};

#endif // FOCUS_INCLUDE_ODEINT_STEPPER_EULER_HPP

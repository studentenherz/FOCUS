#if !defined(FOCUS_INCLUDE_ODEINT_INTEGRATOR_HPP)
#define FOCUS_INCLUDE_ODEINT_INTEGRATOR_HPP

/**
 * Integrate a system with fixed size steps and 
 * method given by the stepper.
 * 
 * @param stepper Stepper that implements the method. Has a do_step method that takes 
 * the arguments `do_step(system, state, current_time, delta_time)`; calculates the state at
 * current_time + delta_time and sets state to it.
 * @param sys Equation system. It's called with arguments `(x, dxdy, t)`, `x` is the current
 * state and `t` the time; sets `dxdy` to the derivative of x with respect to t at current state.
 * @param x State that will be updated during integration.
 * @param t0 Initial time
 * @param dt Time delta
 * @param Nsteps Amount of steps of integration
 * @param obs Observes the state during the integration
 * @param obs_skip_steps Amount of steps to skip between observations (default = 0)
 */
template<typename stepper_type, typename system_type, typename state_type, typename scalar_type, typename observer_type>
#ifdef __CUDACC__
__host__ __device__
#endif
size_t integrate(stepper_type& stepper, system_type& sys, state_type& x, scalar_type t0, scalar_type dt, size_t Nsteps, observer_type& obs, size_t obs_skip_steps = 0){
	scalar_type t = t0;
	size_t step = 0;
	while(step < Nsteps){
		if (step % (obs_skip_steps + 1) == 0)
			obs(x, t);

		stepper.do_step(sys, x, t, dt);
		if(hasnan(x))
			break;
		t += dt;
		step++;
	}

	return step;
}

/**
 * Integrate a system with fixed size steps and 
 * method given by the stepper.
 * 
 * @param stepper Stepper that implements the method. Has a do_step method that takes 
 * the arguments `do_step(system, state, current_time, delta_time)`; calculates the state at
 * current_time + delta_time and sets state to it.
 * @param sys Equation system. It's called with arguments `(x, dxdy, t)`, `x` is the current
 * state and `t` the time; sets `dxdy` to the derivative of x with respect to t at current state.
 * @param x State that will be updated during integration.
 * @param t0 Initial time
 * @param dt Time delta
 * @param Nsteps Amount of steps of integration
 */
template<typename stepper_type, typename system_type, typename state_type, typename scalar_type>
#ifdef __CUDACC__
__host__ __device__
#endif
size_t integrate(stepper_type& stepper, system_type& sys, state_type& x, scalar_type t0, scalar_type dt, size_t Nsteps){
scalar_type t = t0;
	size_t step = 0;
	while(step < Nsteps){
		stepper.do_step(sys, x, t, dt);
		if(hasnan(x))
			break;
		t += dt;
		step++;
	}

	return step;
}

/**
 * Integrate a system with fixed size steps and 
 * method given by the stepper.
 * 
 * @param stepper Stepper that implements the method. Has a do_step method that takes 
 * the arguments `do_step(system, state, current_time, delta_time)`; calculates the state at
 * current_time + delta_time and sets state to it.
 * @param sys Equation system. It's called with arguments `(x, dxdy, t)`, `x` is the current
 * state and `t` the time; sets `dxdy` to the derivative of x with respect to t at current state.
 * @param x State that will be updated during integration.
 * @param t0 Initial time
 * @param dt Time delta
 * @param Nsteps Amount of steps of integration
 */
template<typename stepper_type, typename system_type, typename state_type, typename scalar_type, typename stopping_condition_t>
#ifdef __CUDACC__
__host__ __device__
#endif
size_t stopping_condition_integrate(stepper_type& stepper, system_type& sys, state_type& x, scalar_type t0, scalar_type dt, size_t Nsteps, stopping_condition_t stopping_condition){
scalar_type t = t0;
	size_t step = 0;
	while(step < Nsteps){
		if(hasnan(x) || stopping_condition(sys, x, t))
			break;
		stepper.do_step(sys, x, t, dt);
		t += dt;
		step++;
	}

	return step;
}

#endif // FOCUS_INCLUDE_ODEINT_INTEGRATOR_HPP

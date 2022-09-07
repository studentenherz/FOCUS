#if !defined(FOCUS_INCLUDE_ODEINT_STEPPER_COLLISION_STEPPER_HPP)
#define FOCUS_INCLUDE_ODEINT_STEPPER_COLLISION_STEPPER_HPP

template<typename system_type, typename state_type, typename scalar_type, typename orbit_stepper_type, typename collision_operator_type>
class CollisionStepper{
	const size_t collisions_nstep;
	size_t steps = 0;
	orbit_stepper_type orbit_stepper;
	collision_operator_type collision_operator;
public:
	CollisionStepper(size_t nstep, collision_operator_type collisions): collisions_nstep(nstep), collision_operator(collisions) {}
	
	void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		orbit_stepper.do_step(sys, x, t, dt);
		if (++steps == collisions_nstep){
			collision_operator.euler_step(x, t, dt * collisions_nstep);
			steps = 0;
		}
	}
};

#endif // FOCUS_INCLUDE_ODEINT_STEPPER_COLLISION_STEPPER_HPP

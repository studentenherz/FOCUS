#if !defined(FOCUS_INCLUDE_ODEINT_STEPPER_RK46_NL_HPP)
#define FOCUS_INCLUDE_ODEINT_STEPPER_RK46_NL_HPP

/**
 * [Low-dissipation and low-dispersion fourth-order Runge–Kutta](https://acoustique.ec-lyon.fr/publi/berland_cfluids06.pdf)
 */
template<typename system_type, typename state_type, typename scalar_type>
class RK46NL{
	scalar_type a[6] = {0.000000000000, -0.737101392796, -1.634740794341, -0.744739003780, -1.469897351522, -2.813971388035};
	scalar_type b[6] = {0.032918605146,  0.823256998200,  0.381530948900,  0.200092213184,  1.718581042715,  0.270000000000};
	scalar_type c[6] = {0.000000000000,  0.032918605146,  0.249351723343,  0.466911705055,  0.582030414044,  0.847252983783};
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	RK46NL() {}
	
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		state_type dxdt, xx; // this assumes default constructor gives a zero state

		for (size_t i = 0; i < 6; i++){
			scalar_type tt = t + c[i] * dt;
			sys(x, dxdt, tt); // dxdt = f(x, t)

			xx = a[i] * xx + dxdt * dt;
			x = x + b[i] * xx;
		}

	}
};

#endif // FOCUS_INCLUDE_ODEINT_STEPPER_RK46_NL_HPP

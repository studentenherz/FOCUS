#include <fstream>
#include <iostream>
#include <string>

#include "collisions/focker_plank.hpp"
#include "lorentz.hpp"
#include "odeint/integrator.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "types/equilibrium.hpp"
#include "magnetic_field.hpp"
#include "formats/geqdsk.hpp"
#include "formats/input_gacode.hpp"
#include "cxxopts.hpp"


template<typename system_type, typename state_type, typename scalar_type, typename collision_op_t>
class CollisionStepper{
	const size_t collisions_nstep;
	size_t steps = 0;
	RK46NL<system_type, state_type, scalar_type> orbit_stepper;
	collision_op_t& collision_operator;
public:
	CollisionStepper(size_t nstep, collision_op_t& collisions): collisions_nstep(nstep), collision_operator(collisions) {}
	
	void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		orbit_stepper.do_step(sys, x, t, dt);
		if (++steps == collisions_nstep){
			collision_operator.euler_step(x, t, dt * collisions_nstep);
			steps = 0;
		}
	}
};

class FileObserver{
	std::ofstream &_fo;
	double _a, _v0, _Omega;
	bool cart;
public:
	FileObserver(std::ofstream& fo, double a, double v0, double Omega,bool cart_coord = false): _fo(fo), _a(a), _v0(v0), _Omega(Omega), cart(cart_coord) {}
	
	void operator()(State v, double t){
		if (cart)
			_fo << t / _Omega << '\t' << v[0] * cos(v[1]) * _a << ' ' << v[0] * sin(v[1]) * _a << ' ' << v[2] * _a << ' ' << (v[3] * cos(v[1]) - v[4] * sin(v[1])) * _v0 << ' ' << (v[3] * sin(v[1]) + v[3] * cos(v[1])) * _v0 << ' ' << v[5] * _v0<< '\n';
		else
			_fo << t/ _Omega << '\t' << v[0]  * _a << ' ' << v[1]  << ' ' << v[2] * _a << ' ' << v[3] * _v0 << ' ' << v[4] * _v0 << ' ' << v[5] * _v0 << '\n';
	}
};

/**
 * A vector field class that returns zero
 */
class HomogeneusNullB {
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	HomogeneusNullB(){}
	
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Vector3 operator()(Vector3 /* r */, double /* t */) {Vector3 zero; return zero;}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double psi(Vector3 /* r */, double /* t */ ){
		return 0.05;
	}
};

int main(int argc, char* argv[]){
	cxxopts::options options(argv[0], "Simulate particle with only slowing down collisions");

	options.add_options()
		("i,input", "input.gacode input file", cxxopts::value<std::string>())
		("h,help", "Show this help message")
	;

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			std::cout << options.help() << "\n";
			return EXIT_SUCCESS;
		}

		double pre_eta = 2.42569424e18;
		double kappa = 0.0168686;
		double logl_prefactor = 18.4527;
		
		double q_over_m =  9.64853e7; // C/kg  (e/1 Da)
		double Omega = q_over_m; // cyclotron frequency
		double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
		double a = 2; // m
		double gam = v0 / (a * Omega); // dimensionless factor

		double eta = pre_eta / (Omega * v0 * v0 * v0);
		std::cout << "eta: " << eta << '\n';


		// Particles
		Particle alpha(2.01410177811, 1); // Deuteron

		std::vector<std::string> species_identifiers;
		Plasma plasma = read_input_gacode(result["input"].as<std::string>(), species_identifiers);
		plasma.logl_prefactor = logl_prefactor;


		HomogeneusNullB B;

		std::cout << plasma[0].logl(alpha, 0.05, 0) << '\n';

		// System with Lorentz force
		typedef Lorentz<NullForce, HomogeneusNullB, NullVectorField> System;
		System sys(gam, alpha, B, null_vector_field, null_force);

		double dt = 0.1;
		size_t N = 810000000;
		size_t n_obs =   999;

		for (unsigned long long seed = 1; seed < 2; seed++){

			NormalRand ran(seed);

			// Collisions operator
			typedef FockerPlank<NormalRand, HomogeneusNullB> collision_t;
			collision_t collisions(plasma, alpha, B, eta, kappa, ran);

			// Stepper
			CollisionStepper<System, State, double, collision_t> stepper(200, collisions);

			// Initial step
			State x = {2.0 / a, 0, 0, 0.0, 0.0, 0.1621};

			std::string fname = "slow_down/" + std::to_string(seed) + ".dat";
			// Observer
			std::ofstream fo(fname);
			FileObserver obs(fo, a, v0, Omega, true);

			std::cout << "Calculating " << fname << '\n';
			integrate(stepper, sys, x, 0.0, dt, N, obs, n_obs);

			fo.close();
		}

		return EXIT_SUCCESS;
	}
	catch(cxxopts::option_error const& e){
		std::cerr << e.what() << '\n';
		return EXIT_FAILURE;
	}
}
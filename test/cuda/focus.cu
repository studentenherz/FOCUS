#include <iostream>
#include <fstream>
#include <cstdlib>

#include "cuda/paralell_reduction.hpp"
#include "formats/particle_states.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "odeint/stepper/boris.hpp"
#include "formats/input_gacode.hpp"
#include "types/equilibrium.hpp"
#include "odeint/integrator.hpp"
#include "formats/geqdsk.hpp"
#include "magnetic_field.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "types/plasma.hpp"
#include "collisions/focker_plank.hpp"
#include "lorentz.hpp"
#include "cxxopts.hpp"
#include "formats/matrix.hpp"

void adimensionalize(Array<State> states, double v0, double a){
	for (size_t i = 0; i < states.size(); i++){
		states[i][0] /= a;
		states[i][2] /= a;
		states[i][3] /= v0;
		states[i][4] /= v0;
		states[i][5] /= v0;
	}
	// dump_states("dimensionless_states.dat", states, Particle(1, 1));
}

void dimensionalize(Array<State> states, double v0, double a){
	for (size_t i = 0; i < states.size(); i++){
		states[i][0] *= a;
		states[i][2] *= a;
		states[i][3] *= v0;
		states[i][4] *= v0;
		states[i][5] *= v0;
	}
}

template<typename system_type, typename state_type, typename scalar_type, typename CollisionOperator_t>
class CollisionStepper{
	const size_t collisions_nstep;
	size_t steps;
	// Boris<system_type, state_type, scalar_type> orbit_stepper;
	RK46NL<system_type, state_type, scalar_type> orbit_stepper;
	CollisionOperator_t& collision_operator;
public:
	__host__ __device__
	CollisionStepper(size_t offset, size_t nstep, CollisionOperator_t& collisions): collisions_nstep(nstep), steps(offset % collisions_nstep), collision_operator(collisions) {}

	__host__ __device__	
	void do_step(system_type sys, state_type& x, scalar_type t, scalar_type dt){
		orbit_stepper.do_step(sys, x, t, dt);
		if (++steps == collisions_nstep){
			collision_operator.euler_step(x, t, dt * collisions_nstep);
			steps = 0;
		}
	}
};

typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> system_t;

__global__
void k_integrate(MagneticFieldMatrix B_matrix, Equilibrium eq, Plasma plasma, Particle test_particle, double gamma, double eta, double kappa, PhiloxCuRand philox, Array<State> x, double t0, double dt, size_t Nsteps, size_t offset, Array<double> pitch, Array<double> energy){
	// Index of thread	
	size_t idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < x.size()){
		// Magnetic field
		MagneticFieldFromMatrix B(B_matrix, eq.bcentr);
		
		// Lorentz Force Equiations System
		system_t sys(gamma, test_particle, B, d_null_vector_field, d_null_force);
		
		// Collision operator
		FockerPlank<PhiloxCuRand, MagneticFieldFromMatrix> collision(plasma, test_particle, B, eta, kappa, philox);

		// Stepper
		CollisionStepper<system_t, State, double, FockerPlank<PhiloxCuRand, MagneticFieldFromMatrix>> stepper(offset, 200, collision);

		// Integrate
		integrate(stepper, sys, x[idx], t0, dt, Nsteps);

		Vector3 r = get_position(x[idx]);
		Vector3 v = get_velocity(x[idx]);
		Vector3 b = B(r, 0);

		pitch[idx] = pitch_between(v, b);
		energy[idx] = dot(v, v);
	}
}

bool make_dir(std::string dirname){
	std::string command = "mkdir -p ";
	command += dirname;
	return  (system(command.c_str()) == 0);
}

void device_integrate(MagneticFieldMatrix& B_matrix, Equilibrium& eq, Plasma& plasma, Particle test_particle, double gamma, double eta, double kappa, Array<State>& states, double t0, double dt, size_t Nsteps, size_t Nobs, std::string out_dir, double v0, double a){
	// Random generator
	PhiloxCuRand philox(states.size());

	static const size_t blockSize = 1024;
	static const size_t gridSize = 24;

	kernel_init_philox_rand<<<gridSize, blockSize>>>(philox, 1234ull);

	if (make_dir(out_dir))
		std::cout << "\x1b[32mSuccess\x1b[0m creating output directory \x1b[1m" << out_dir << "\x1b[0m\n";
	else{
		std::cerr << "\x1b[31mError\x1b[0m creating output directory \x1b[1m" << out_dir << "\x1b[0m!\n";
		exit(-1);
	}

	std::ofstream index_file(out_dir + "/.focus_index");
	if (!index_file.is_open()){
		std::cerr << "\x1b[31mError\x1b[0m creating index file\n";
		exit(-1);
	}
	std::string curr_file;

	adimensionalize(states, v0, a);
	
	// Store pitch and energy in host
	Array<double> pitch(states.size());
	Array<double> energy(states.size());

	// Same for device
	Array<double> d_pitch; d_pitch.construct_in_host_for_device(pitch);
	Array<double> d_energy; d_energy.construct_in_host_for_device(energy);

	Array<State> d_states;
	d_states.construct_in_host_for_device(states);

	Plasma d_plasma(0, 0, 0);
	d_plasma.construct_in_host_for_device(plasma);

	Equilibrium d_eq;
	d_eq.construct_in_host_for_device(eq);

	MagneticFieldMatrix d_B_matrix;
	d_B_matrix.construct_in_host_for_device(B_matrix);


	double t = t0;
	for (size_t i = 0; i <= Nsteps / Nobs; i++){
		k_integrate<<<gridSize, blockSize>>>(d_B_matrix, d_eq, d_plasma, test_particle, gamma, eta, kappa, philox, d_states, t, dt, Nobs, i * Nobs, d_pitch, d_energy);
		states.copy_to_host_from_device(d_states);
		pitch.copy_to_host_from_device(d_pitch);
		energy.copy_to_host_from_device(d_energy);
		t += Nobs * dt;

		std::cout << "Saving states at " << t;
		dimensionalize(states, v0, a);
		curr_file = "states_" + std::to_string(t) + ".dat";
		
		if (dump_states(out_dir + "/" + curr_file, states, test_particle, "pitch | E[dimensionless]" ,pitch, energy)){
			index_file << curr_file << std::endl;
			std::cout << "\x1b[32m done\x1b[0m\n";
		}
	}

	std::cout << "\x1b[32mAll done\x1b[0m, no errors (that I know of)\n";
}

int main(int argc, char* argv[]){
	cxxopts::options options(argv[0], "Test everything together");

	options.add_options()
		("g,geqdsk", "G-EQDSK file", cxxopts::value<std::string>())
		("i,input-gacode", "input.gacode file", cxxopts::value<std::string>())
		("d,output-dir", "Output directory", cxxopts::value<std::string>()->default_value("focus_output"))
		("s,initial-states", "File with initial states", cxxopts::value<std::string>())
		("h,help", "Show this help message")
	;

	try{
		auto result = options.parse(argc, argv);
	
		if (result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		Equilibrium eq = read_geqdsk(result["geqdsk"].as<std::string>());
		MagneticFieldMatrix B_matrix(eq, 26, 600, false);

		double t0 = 0;
		double dt = 0.1;
		// size_t N 	= 5000000;
		size_t N =   10000; // Enough for slow down
		size_t nobs = 500;
		// size_t N = 410000000; // Enough for slow down
		// size_t nobs = 1000000;
		// size_t N = 10;
		// size_t nobs = 1;
		// size_t N 	= 500000;
		// size_t nobs = 50000;
		double logl_prefactor = 18.4527;

		double q_over_m =  9.64853e7; // C/kg electron charge / Da
		double Omega = q_over_m * eq.bcentr; // cyclotron frequency
		double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
		double a = eq.rdim; // m
		double gamma = v0 / (a * Omega); // dimensionless factor

		std::cout << 1/Omega << '\n';
	
		double eta = 2.27418e-12;
		double kappa = 0.016870543;


		Plasma plasma = read_input_gacode(result["input-gacode"].as<std::string>());
		plasma.logl_prefactor = logl_prefactor;
		Particle test_particle(1, 2.01410177811);


		Array<State> states = load_states(result["initial-states"].as<std::string>());
		device_integrate(B_matrix, eq, plasma, test_particle, gamma, eta, kappa, states, t0, dt, N, nobs, result["output-dir"].as<std::string>(), v0, a);
		return 0;
		
	}
	catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
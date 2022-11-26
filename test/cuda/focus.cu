#include <iostream>
#include <fstream>
#include <cstdlib>

#include "collisions/atomic_processes.hpp"
#include "collisions/focker_plank.hpp"
#include "cuda/paralell_reduction.hpp"
#include "formats/particle_states.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "odeint/stepper/boris.hpp"
#include "formats/input_gacode.hpp"
#include "handle_cuda_errors.hpp"
#include "types/equilibrium.hpp"
#include "odeint/integrator.hpp"
#include "magnetic_field.hpp"
#include "formats/geqdsk.hpp"
#include "formats/matrix.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "types/plasma.hpp"
#include "lorentz.hpp"
#include "cxxopts.hpp"

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

template<typename system_type, typename state_type, typename scalar_type, typename CollisionOperator_t, typename RandomGenerator_t, typename MagneticField_t>
class CollisionStepper{
	const size_t collisions_nstep;
	const size_t atomic_nstep;
	size_t steps;
	Boris<system_type, state_type, scalar_type> boris;
	RK46NL<system_type, state_type, scalar_type> rk46nl;
	CollisionOperator_t& collision_operator;
	// AtomicProcesses_t& atomic_processes;

	Array<AtomicProcess>& processes;
	RandomGenerator_t& ran_gen;
	MagneticField_t& B;
	Plasma& plasma;
	double energy_conversion_factor;
	double n_conversion_factor_4_APs;
	double Omega;
public:
	__host__ __device__
	CollisionStepper(size_t offset, size_t nstep, CollisionOperator_t& collisions, size_t atomic_nstep, Array<AtomicProcess>& processes, RandomGenerator_t& ran_gen, MagneticField_t& B, Plasma& plasma, double energy_conversion, double n_conversion, double Omega): collisions_nstep(nstep), steps(offset), collision_operator(collisions), atomic_nstep(atomic_nstep), processes(processes), ran_gen(ran_gen), B(B), plasma(plasma), energy_conversion_factor(energy_conversion), n_conversion_factor_4_APs(n_conversion), Omega(Omega) {}

	__host__ __device__	
	void do_step(system_type& sys, state_type& x, scalar_type t, scalar_type dt){
		steps++;
		
		/*
			Boris stepper is very good for particles subjected to Lorentz force,
			but for neutral particles is not very good. Use RK46NL instead in that case.
		*/
		if (sys.part.q == 0){
			rk46nl.do_step(sys, x, t, dt);

			
			dt /= Omega; // Time in seconds
			sys.part.t += dt;

			Vector3 r = get_position(x);
			Vector3 v = get_velocity(x);

			double energy = sys.part.m * dot(v, v) * energy_conversion_factor;
			double psi = B.psi(r, t);

			Array<double> p(processes.size());
			for (size_t i = 0; i < processes.size(); i++)
				p[i] = processes[i].P(sys.part, energy, psi, plasma, t, dt, n_conversion_factor_4_APs);
			
			for (size_t i = 1; i < processes.size(); i++)
				p[i] += p[i - 1];

			// if (threadIdx.x + blockDim.x * blockIdx.x == 0)
			// 	printf("P = %e\n", p[p.size() - 1]);

			double ran = ran_gen.uniform();

			for (size_t i = 0 ; i < processes.size(); i++)
				if (ran < p[i]){
					processes[i].apply(sys.part);
					return;
				}

			// if (steps % atomic_nstep == 0)
			// 	atomic_processes(sys.part, x, t, dt * atomic_nstep);
		}
		else{
			boris.do_step(sys, x, t, dt);

			if (steps % collisions_nstep == 0)
				collision_operator.euler_step(x, t, dt * collisions_nstep);
		}
	}
};

template<typename System_t>
class IonizedStoppingCondition{
	double rbdry_min, rbdry_max;
	double zbdry_min, zbdry_max;
	State& birth;
public:
	__device__
	IonizedStoppingCondition(double rbdry_min, double rbdry_max, double zbdry_min, double zbdry_max, State& birth): rbdry_min(rbdry_min), rbdry_max(rbdry_max), zbdry_min(zbdry_min), zbdry_max(zbdry_max), birth(birth) {}
	__device__
	bool operator()(System_t& sys, State& x, double t){
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (sys.part.q != 0){
			birth = x;
			printf("Ionized particle at thread %d\n", idx);
			sys.part.n = 777;
			return true;
		}
		// Check if it's out and not comig in
		if (x[0] < rbdry_min || (x[0] > rbdry_max && x[3] > 0) || (x[2] < zbdry_min && x[5] < 0) || (x[2] > zbdry_max && x[5] > 0)){
			printf("Out of boundary (%lf, %lf, %lf) particle at thread %d\n",x[0], x[1], x[2], idx);
			sys.part.n = 777;
			return true;
		}
		return false;
	}
};

typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> system_t;

__global__
void k_integrate(MagneticFieldMatrix B_matrix, Equilibrium eq, Plasma plasma, Array<Particle> test_particles, double gamma, double eta, double kappa, PhiloxCuRand philox, Array<State> x, double t0, double dt, size_t Nsteps, size_t offset, Array<double> pitch, Array<double> energy, Array<AtomicProcess> atomic_processes, double Omega, double energy_conversion_factor, double n_conversion_factor_4_APs, State initial_injection_state, double rbdry_min, double rbdry_max, double zbdry_min, double zbdry_max, Array<State> birth){
	// Index of thread	
	size_t idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < x.size()){
		if (test_particles[idx].n == 777) return;

		// Magnetic field
		MagneticFieldFromMatrix B(B_matrix, eq.bcentr);
		
		// Lorentz Force Equiations System
		system_t sys(gamma, test_particles[idx], B, d_null_vector_field, d_null_force);
		
		// Collision operator
		FockerPlank<PhiloxCuRand, MagneticFieldFromMatrix> collision(plasma, test_particles[idx], B, eta, kappa, philox);

		// AtomicProcessesHandler<PhiloxCuRand, MagneticFieldFromMatrix> atomic_processes_handler(atomic_processes, philox, B, plasma, Omega, energy_conversion_factor, n_conversion_factor_4_APs);

		// Stepper
		CollisionStepper<system_t, State, double, FockerPlank<PhiloxCuRand, MagneticFieldFromMatrix>, PhiloxCuRand, MagneticFieldFromMatrix> stepper(offset, 200, collision, 1, atomic_processes, philox, B, plasma, energy_conversion_factor, n_conversion_factor_4_APs, Omega);

		IonizedStoppingCondition<system_t> isc(rbdry_min, rbdry_max, zbdry_min, zbdry_max, birth[idx]);

		// Integrate
		stopping_condition_integrate(stepper, sys, x[idx], t0, dt, Nsteps, isc);

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

void device_integrate(MagneticFieldMatrix& B_matrix, Equilibrium& eq, Plasma& plasma, Particle test_particle, double gamma, double eta, double kappa, Array<State>& states, double t0, double dt, size_t Nsteps, size_t Nobs, std::string out_dir, double v0, double a, Array<AtomicProcess>& atomic_processes, double Omega, double energy_conversion_factor, double n_conversion_factor_4_APs){
	// Random generator
	PhiloxCuRand philox(states.size());

	static const size_t blockSize = 1024;
	static const size_t gridSize = 24;

	kernel_init_philox_rand<<<gridSize, blockSize>>>(philox, 1234ull);
	checkKernelErr();

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
	Array<double> d_pitch; gpuErrchk( d_pitch.construct_in_host_for_device(pitch) );
	Array<double> d_energy; gpuErrchk( d_energy.construct_in_host_for_device(energy) );

	Array<State> d_states;
	gpuErrchk( d_states.construct_in_host_for_device(states) );

	Plasma d_plasma(0, 0, 0);
	gpuErrchk( d_plasma.construct_in_host_for_device(plasma) );

	Equilibrium d_eq;
	gpuErrchk( d_eq.construct_in_host_for_device(eq) );

	MagneticFieldMatrix d_B_matrix;
	gpuErrchk( d_B_matrix.construct_in_host_for_device(B_matrix) );

	/*
	 Each of the elements of the array needs to be passed to the device first
	 and then the arrray will be copied with the values of the new pointers to
	 the already in device element.
	*/
	Array<AtomicProcess> pre_device_atomic_processes(atomic_processes.size());
	for (size_t i = 0; i < atomic_processes.size(); i++){
		gpuErrchk( pre_device_atomic_processes[i].construct_in_host_for_device(atomic_processes[i]) );
	}

	Array<AtomicProcess> d_atomic_processes;
	gpuErrchk( d_atomic_processes.construct_in_host_for_device(pre_device_atomic_processes) );

	// Birth
	Array<State> birth(states.size()); // all initialized to 0
	Array<State> d_birth; gpuErrchk( d_birth.construct_in_host_for_device(birth) );

	double rbdry_min = eq.rbdry[0], rbdry_max = eq.rbdry[0];
	double zbdry_min = eq.zbdry[0], zbdry_max = eq.zbdry[0];

	for (size_t i = 1; i < eq.nbdry; i++){
		if (eq.rbdry[i] < rbdry_min) rbdry_min = eq.rbdry[i];
		if (eq.rbdry[i] > rbdry_max) rbdry_max = eq.rbdry[i];
		if (eq.zbdry[i] < zbdry_min) zbdry_min = eq.zbdry[i];
		if (eq.zbdry[i] > zbdry_max) zbdry_max = eq.zbdry[i];
	}


	// Dimensionless is what we need to compare
	rbdry_min /= a;
	rbdry_max /= a;
	zbdry_min /= a;
	zbdry_max /= a;

	std::cout << "r = (" << rbdry_min << ", " << rbdry_max << ")\tz = (" << zbdry_min << ", " << zbdry_max << ")\n";

	State intial_injection_state = states[0];

	Array<Particle> test_particles(states.size());
	for (size_t i = 0; i < test_particles.size(); i++)
		test_particles[i] = test_particle;
	Array<Particle> d_test_particles;
	gpuErrchk( d_test_particles.construct_in_host_for_device(test_particles) );

	double t = t0;
	for (size_t i = 0; i <= Nsteps / Nobs; i++){
		k_integrate<<<gridSize, blockSize>>>(d_B_matrix, d_eq, d_plasma, d_test_particles, gamma, eta, kappa, philox, d_states, t, dt, Nobs, i * Nobs, d_pitch, d_energy, d_atomic_processes, Omega, energy_conversion_factor, n_conversion_factor_4_APs, intial_injection_state, rbdry_min, rbdry_max, zbdry_min, zbdry_max, d_birth);
		cudaDeviceSynchronize();
		checkKernelErr();

		gpuErrchk( states.copy_to_host_from_device(d_states) );
		gpuErrchk( pitch.copy_to_host_from_device(d_pitch) );
		gpuErrchk( energy.copy_to_host_from_device(d_energy) );
		t += Nobs * dt;

		std::cout << "Saving states at " << t;
		dimensionalize(states, v0, a);
		curr_file = "states_" + std::to_string(t) + ".dat";
		
		if (dump_states(out_dir + "/" + curr_file, states, test_particle, "pitch | E[dimensionless]" ,pitch, energy)){
			index_file << curr_file << std::endl;
			std::cout << "\x1b[32m done\x1b[0m\n";
		}
	}

	gpuErrchk( birth.copy_to_host_from_device(d_birth) );
	dimensionalize(birth, v0, a);
	dump_states("birth.dat", birth, test_particle);

	std::cout << "\x1b[32mAll done\x1b[0m, no errors (that I know of)\n";
}

int main(int argc, char* argv[]){
	cxxopts::options options(argv[0], "Test everything together");

	options.add_options()
		("g,geqdsk", "G-EQDSK file", cxxopts::value<std::string>())
		("i,input-gacode", "input.gacode file", cxxopts::value<std::string>())
		("d,output-dir", "Output directory", cxxopts::value<std::string>()->default_value("focus_output"))
		("s,initial-states", "File with initial states", cxxopts::value<std::string>())
		("a,atomic", "Directory with atomic processes", cxxopts::value<std::string>())
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
		double dt = 0.01;
		size_t N 	=   500000;
		size_t nobs =  10000;
		// size_t nobs = 500;
		// size_t N =   10000; // Enough for slow down
		// size_t nobs = 10;
		// size_t N = 410000000; // Enough for slow down
		// size_t nobs = 1000000;
		// size_t N = 1;
		// size_t nobs = 1;
		// size_t N 	= 500000;
		// size_t nobs = 50000;
		double logl_prefactor = 18.4527;

		double q_over_m =  9.64853e7; // C/kg electron charge / Da
		double Omega = q_over_m * eq.bcentr; // cyclotron frequency
		double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
		double a = eq.rdim; // m
		double gamma = v0 / (a * Omega); // dimensionless factor
		

		double energy_conversion_factor = 1.7572e3; // 1/2 m0 * (v0)^2 => keV
		double n_conversion_factor_4_APs = 0.1; // n here is in 10^19 m^-3, atomic processes where calculated for n in 10^20
		double eta = 2.27418e-12;
		double kappa = 0.016870543;

		std::vector<std::string> species_identifiers;
		Plasma plasma = read_input_gacode(result["input-gacode"].as<std::string>(), species_identifiers);
		plasma.logl_prefactor = logl_prefactor;
		Particle test_particle(2.01410177811, 0);

		Array<AtomicProcess> atomic_processes = load_atomic_processes(species_identifiers, result["atomic"].as<std::string>());
		std::cout << "Atomic processes loaded\n";

		Array<State> states = load_states(result["initial-states"].as<std::string>());
		device_integrate(B_matrix, eq, plasma, test_particle, gamma, eta, kappa, states, t0, dt, N, nobs, result["output-dir"].as<std::string>(), v0, a, atomic_processes, Omega, energy_conversion_factor, n_conversion_factor_4_APs);
		return 0;
		
	}
	catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
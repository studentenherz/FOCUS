#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "files.hpp"
#include "formats/geqdsk.hpp"
#include "lorentz.hpp"
#include "magnetic_field.hpp"
#include "types/equilibrium.hpp"
#include "types/vector.hpp"
#include "odeint/integrator.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "cxxopts.hpp"
#include "types/particle.hpp"

class ArrayObserver{
	Array<double> &_times;
	Array<State> &_states;
	double _a, _v0, _Omega;
	bool cart;
	size_t index = 0;
public:
	__host__ __device__
	ArrayObserver(Array<double> &times, Array<State> &states, double a, double v0, double Omega, bool cart_coord = false): _times(times), _states(states), _a(a), _v0(v0), _Omega(Omega), cart(cart_coord) {}
	
	__host__ __device__
	void operator()(State v, double t){
		_times[index] = t / _Omega;
		if (cart)
			_states[index++] = {v[0] * cos(v[1]) * _a, v[0] * sin(v[1]) * _a, v[2] * _a, (v[3] * cos(v[1]) - v[4] * sin(v[1])) * _v0, (v[3] * sin(v[1]) + v[3] * cos(v[1])) * _v0, v[5] * _v0};
		else
			_states[index++] = {v[0]  * _a, v[1] , v[2] * _a, v[3] * _v0, v[4] * _v0, v[5] * _v0};
	}
};

// Integration Kernel
__global__
void k_integrate(MagneticFieldMatrix B_matrix, Equilibrium eq, State x0, double t0, double dt, size_t Nsteps, Array<double> times, Array<State> states, size_t nskip){
	MagneticFieldFromMatrix B(B_matrix, eq.bcentr);

	double q_over_m =  9.58e7; // C/kg proton
	double Omega = q_over_m * eq.bcentr; // cyclotron frequency
	double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
	double a = eq.rdim; // m
	double gam = v0 / (a * Omega); // dimensionless factor

	Particle part(1, 1.04);
	typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> System;
	System sys(gam, part, B, d_null_vector_field, d_null_force);
	State x = x0;
	RK46NL<System, State, double> rk46nl;
	ArrayObserver obs(times, states, a, v0, Omega, true);

	integrate(rk46nl, sys, x, t0, dt, Nsteps, obs, nskip);
}

void integrate_in_device(MagneticFieldMatrix& B_matrix, Equilibrium& eq, State x0, double t0, double dt, size_t Nsteps, std::string ofname, size_t nskip){
	size_t Nout = size_t(Nsteps / (nskip + 1));
	Array<State> h_states(Nout);
	Array<State> d_states;
	d_states.construct_in_host_for_device(h_states);

	Array<double> h_times(Nout);
	Array<double> d_times;
	d_times.construct_in_host_for_device(h_times);

	k_integrate<<<1, 1>>>(B_matrix, eq, x0, t0, dt, Nsteps, d_times, d_states, nskip);

	h_states.copy_to_host_from_device(d_states);
	h_times.copy_to_host_from_device(d_times);

	std::ofstream fo(ofname);

	for(size_t i =0; i < h_states.size(); i++){
		fo << h_times[i] << ' ' << h_states[i] << '\n';
	}

	fo.close();
}

void integrate_in_host(MagneticFieldMatrix& B_matrix, Equilibrium& eq, State x0, double t0, double dt, size_t Nsteps, std::string ofname, size_t nskip){
	MagneticFieldFromMatrix B(B_matrix, eq.bcentr);

	double q_over_m =  9.58e7; // C/kg proton
	double Omega = q_over_m * eq.bcentr; // cyclotron frequency
	double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
	double a = eq.rdim; // m
	double gam = v0 / (a * Omega); // dimensionless factor

	Particle part(1, 1.04);
	typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> System;
	System sys(gam, part, B, null_vector_field, null_force);
	State x = x0;
	RK46NL<System, State, double> rk46nl;

	size_t Nout = size_t(Nsteps / (nskip + 1));
	Array<State> states(Nout);
	Array<double> times(Nout);

	ArrayObserver obs(times, states, a, v0, Omega, true);
	integrate(rk46nl, sys, x, t0, dt, Nsteps, obs, nskip);

	std::ofstream fo(ofname);

	for(size_t i =0; i < states.size(); i++){
		fo << times[i] << ' ' << states[i] << '\n';
	}

	fo.close();
}

void output_b(MagneticFieldMatrix& B_matrix, Equilibrium& eq){
	dump("Br.dat", B_matrix.Br, false);
	dump("Bt.dat", B_matrix.Bt, false);
	dump("Bz.dat", B_matrix.Bz, false);


	std::cout << B_matrix.r_min * eq.rdim << ' ' << B_matrix.r_max * eq.rdim << ' ' << B_matrix.z_min * eq.rdim << ' ' << B_matrix.z_max * eq.rdim << '\n';
	std::ofstream fobdry("bdry.dat");

	for (size_t i = 0; i < eq.nbdry; i++)
		fobdry << eq.rbdry[i] << ' ' << eq.zbdry[i] << '\n';
	fobdry.close();

	std::ofstream folim("lim.dat");
	for (size_t i = 0; i < eq.nlim; i++)
		folim << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';
	folim.close();
}

int main(int argc, char* argv[]){
	cxxopts::options options("integrate", "Test integrate in CUDA");

	options.add_options()
		("input_file", "", cxxopts::value<std::string>())
		("host_file", "", cxxopts::value<std::string>()->default_value("host.dat"))
		("device_file", "", cxxopts::value<std::string>()->default_value("device.dat"))
		("b,magnetic-field", "Output magnetic field data")
		("h,help", "Show this help message");

	options.positional_help("<G-EQDSK input file> [<Host integration out file> [<Device integration out file>]]");
	options.parse_positional({"input_file", "host_file", "device_file"});

	try{
		auto result = options.parse(argc, argv);
	
		if (result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		Equilibrium eq = read_geqdsk(result["input_file"].as<std::string>());
		MagneticFieldMatrix B_matrix(eq, 26, 600);

		if (result.count("magnetic-field")){
			output_b(B_matrix, eq);
			return 0;
		}

		State x0 = {2.2 / eq.rdim, 0, 0, 0.0, 0.01, 0.3};
		double t0 = 0;
		double dt = 0.0001;
		size_t N 	= 3100000;
		size_t nskip 	= 999;
		
		integrate_in_host(B_matrix, eq, x0, t0, dt, N, result["host_file"].as<std::string>(), nskip);
		integrate_in_device(B_matrix, eq, x0, t0, dt, N, result["device_file"].as<std::string>(), nskip);

		return 0;
	}
	catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
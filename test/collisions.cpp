#include <fstream>
#include <iostream>
#include <string>

#include "collisions.hpp"
#include "lorentz.hpp"
#include "odeint/integrator.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/vector.hpp"
#include "types/equilibrium.hpp"
#include "magnetic_field.hpp"
#include "geqdsk.hpp"


template<typename system_type, typename state_type, typename scalar_type>
class CollisionStepper{
	const size_t collisions_nstep;
	size_t steps = 0;
	RK46NL<system_type, state_type, scalar_type> orbit_stepper;
	FockerPlank& collision_operator;
public:
	CollisionStepper(size_t nstep, FockerPlank& collisions): collisions_nstep(nstep), collision_operator(collisions) {}
	
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

int main(int argc, char* argv[]){
	if (argc < 2){
		std::cout << "usage:\ncollisions <g-eqdsk file>\n";
		return -1;
	}

	Equilibrium eq = read_geqdsk(argv[1]);
	MagneticFieldMatrix B_matrix(eq, 26, 600);

	MagneticFieldFromMatrix B(B_matrix, eq.bcentr);

	double q_e = 1.0;
	double m_e = 1.0;
	double logl_e = 17.5;
	double eta = 3453.5;
	
	double q_over_m =  9.58e7; // C/kg proton
	double Omega = q_over_m * eq.bcentr; // cyclotron frequency
	double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
	double a = eq.rdim; // m
	double gam = v0 / (a * Omega); // dimensionless factor

	double Tf = 1.61029;
	double nf = 1.0;

	// Particles
	ConstProfileParticle electron(q_e, m_e, logl_e, Tf, nf);
	ConstProfileParticle alpha(1.01, 2.0 * 1836, logl_e, Tf, nf);

	// Particles in plasma
	Array<ParticleSpecies*> plasma(1);
	plasma[0] = &electron;

	// System with Lorentz force
	typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> System;
	System sys(gam, B, null_vector_field, null_force);

	for (unsigned long long seed = 1; seed < 5; seed++){

		// Collisions operator
		FockerPlank collisions(seed, plasma, alpha, eta);

		// Stepper
		CollisionStepper<System, State, double> stepper(200, collisions);

		// Initial step
		State x = {2.2 / a, 0, 0, 0.0, 0.01, 0.3};

		std::string fname = "coll/" + std::to_string(seed) + ".dat";
		// Observer
		std::ofstream fo(fname);
		FileObserver obs(fo, a, v0, Omega, true);

		std::cout << "Calculating " << fname << '\n';
		integrate(stepper, sys, x, 0.0, 0.0001, 30000000, obs, 999);

		fo.close();
	}

	return EXIT_SUCCESS;
}
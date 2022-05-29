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

// Homogeneous Magnetic Field in z direction
struct MagneticField{
	double B0;
	MagneticField(double B): B0(B) {}
	Vector3 operator()(Vector3 /* r */, double /* t */ ){
		Vector3 f;
		f[2] = B0;
		return f;
	}
} B(1);

// Constant temperature profile
double Tf(Vector3 /* r */, double /* t */){
	return 1.61029;
}

// Constant density profile
double nf(Vector3 /* r */, double /* t */){
	return 1.0;
}


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
		if (++steps % collisions_nstep == 0)
			collision_operator.euler_step(x, t, dt * collisions_nstep);
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

int main(){
	double Omega = 0.00085; // cyclotron frequency
	double v0 = 1.29477e7; // m/s (3.54 MeV of a proton)
	double a = 0.5; // m
	double gam = 2.408e6;
	double q_e = 1.0;
	double m_e = 1.0;
	double logl_e = 17.5;
	double eta = 3453.5;

	// Particles
	ParticleSpecies electron(q_e, m_e, logl_e, Tf, nf);
	ParticleSpecies alpha(1.0, 2.0 * 1836, logl_e, Tf, nf);

	// Particles in plasma
	Array<ParticleSpecies> plasma(1);
	plasma[0] = electron;

	// System with Lorentz force
	typedef Lorentz<NullForce, MagneticField, NullVectorField> System;
	System sys(gam, B, null_vector_field, null_force);

	for (unsigned long long seed = 1; seed < 50; seed++){

		// Collisions operator
		FockerPlank collisions(seed, plasma, alpha, eta);

		// Stepper
		CollisionStepper<System, State, double> stepper(200, collisions);

		// Initial step
		State x(1.0, 0.0, 0.0, 0.0, 0.0, 0.23057);

		std::string fname = "coll/" + std::to_string(seed) + ".dat";
		// Observer
		std::ofstream fo(fname);
		FileObserver obs(fo, a, v0, Omega, true);

		std::cout << "Calculating " << fname << '\n';
		integrate(stepper, sys, x, 0.0, 0.00001, 300000, obs, 999);

		fo.close();
	}

	return EXIT_SUCCESS;
}
#include <iostream>
#include <fstream>

#include "files.hpp"
#include "geqdsk.hpp"
#include "lorentz.hpp"
#include "magnetic_field.hpp"
#include "types/equilibrium.hpp"
#include "types/vector.hpp"
#include "odeint/integrator.hpp"
#include "odeint/stepper/rk46_nl.hpp"


class FileObserver{
	std::ofstream &_fo;
	double _a, _v0, _Omega;
public:
	FileObserver(std::ofstream& fo, double a, double v0, double Omega): _fo(fo), _a(a), _v0(v0), _Omega(Omega) {}
	
	void operator()(State v, double t){
		_fo << t / _Omega << '\t' << v[0] * cos(v[1]) * _a << ' ' << v[0] * sin(v[1]) * _a << ' ' << v[2] * _a << ' ' << (v[3] * cos(v[1]) - v[4] * sin(v[1])) * _v0 << ' ' << (v[3] * sin(v[1]) + v[3] * cos(v[1])) * _v0 << ' ' << v[5] * _v0<< '\n';
	}
};

int main(int argc, char* argv[]){
	if (argc < 3){
		std::cout << "usage:\nmagnetic_field <g-eqdsk file> <output_file>\n";
		return -1;
	}

	
	Equilibrium eq = read_eqdsk(argv[1]);
	MagneticFieldMatrix B_matrix(eq, 26, 600);
	MagneticField B(B_matrix, eq.bcentr);

	double q_over_m =  -1.75882001076e10; // C/kg electron
	double Omega = q_over_m * eq.bcentr; // cyclotron frequency
	double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
	double a = eq.rdim; // m
	double gam = v0 / (a * Omega); // dimensionless factor

	std::cout << "gam = " << gam << '\n';

	typedef Lorentz<NullForce, MagneticField, NullVectorField> System;
	System sys(gam, B, null_vector_field, null_force);
	State x(2.0 / a, 0.0, 0.0, 0.0, 1.6, 0.0);
	RK46NL<System, State, double> rk46nl;

	std::ofstream fo(argv[2]);
	FileObserver obs(fo, a, v0, Omega);

	// dump("Br.dat", B_matrix.Br, false);
	// dump("Bt.dat", B_matrix.Bt, false);
	// dump("Bz.dat", B_matrix.Bz, false);
	// dump("ch_psi.dat", B_matrix.ch_psi, false);

	integrate(rk46nl, sys, x, 0.0, 0.001, 1000000, obs, 9);

	return 0;
}
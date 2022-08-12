#include <iostream>
#include <fstream>
#include <stdlib.h>

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
	double _a, _v0, _Omega, _B0;
	bool cart;
	MagneticFieldFromMatrix& _B;
	FineEquilibrium& fineq;
public:
	FileObserver(std::ofstream& fo, double a, double v0, double Omega, MagneticFieldFromMatrix& B, FineEquilibrium& eq, bool cart_coord = false): _fo(fo), _a(a), _v0(v0), _Omega(Omega),  _B0(B.B0()), cart(cart_coord), _B(B), fineq(eq) {}
	
	void operator()(State v, double t){
		// magnetic field 
		Vector3 B = _B(get_position(v), t);
		if (cart)
			_fo << t / _Omega << '\t' << v[0] * cos(v[1]) * _a << ' ' << v[0] * sin(v[1]) * _a << ' ' << v[2] * _a << ' ' << (v[3] * cos(v[1]) - v[4] * sin(v[1])) * _v0 << ' ' << (v[3] * sin(v[1]) + v[3] * cos(v[1])) * _v0 << ' ' << v[5] * _v0 << ' ' << (B[0] * cos(v[1]) - B[1] * sin(v[1])) * _B0 << ' ' << (B[0] * sin(v[1]) + B[1] * cos(v[1])) * _B0 << ' ' << B[2] * _B0 << ' ' << fineq.Psi(v[0], v[2]) << ' ' << fineq.F(v[0], v[2]) << '\n';
		else
			_fo << t/ _Omega << '\t' << v[0]  * _a << ' ' << v[1]  << ' ' << v[2] * _a << ' ' << v[3] * _v0 << ' ' << v[4] * _v0 << ' ' << v[5] * _v0 << ' ' << B * _B0 << ' ' << fineq.Psi(v[0], v[2]) << ' ' << fineq.F(v[0], v[2]) << '\n';
	}
};

int main(int argc, char* argv[]){
	if (argc < 3){
		std::cout << "usage:\nmagnetic_field <g-eqdsk file> <output_file>\n";
		return -1;
	}

	
	Equilibrium eq = read_geqdsk(argv[1]);
	FineEquilibrium fineq(eq, 26);
	MagneticFieldMatrix B_matrix(eq, 26, 600);
	
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
	
	std::cout << "Done\n";

	MagneticFieldFromMatrix B(B_matrix, eq.bcentr);
	// MagneticField B(fineq);

	double q_over_m =  9.58e7; // C/kg proton
	double Omega = q_over_m * eq.bcentr; // cyclotron frequency
	double v0 = 1.84142e7; // m/s (3.54 MeV of a proton)
	double a = eq.rdim; // m
	double gam = v0 / (a * Omega); // dimensionless factor

	typedef Lorentz<NullForce, MagneticFieldFromMatrix, NullVectorField> System;
	System sys(gam, B, null_vector_field, null_force);
	State x = {2.2 / a, 0, 0, 0.0, 0.01, 0.3};
	RK46NL<System, State, double> rk46nl;

	std::ofstream fo(argv[2]);


	FileObserver obs(fo, a, v0, Omega, B, fineq, false);
	integrate(rk46nl, sys, x, 0.0, 0.0001, 30000000, obs, 999);

	return 0;
}
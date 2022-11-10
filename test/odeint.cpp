#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include "odeint/integrator.hpp"
#include "odeint/stepper/euler.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "odeint/stepper/boris.hpp"
#include "types/vector.hpp"
#include "lorentz.hpp"

double gam = 0.0;

// struct System
// {
// 	System() {}
// 	void operator()(const Vector3& x, Vector3& dxdt, double t){
// 		dxdt[0] = x[1];
// 		dxdt[1] = -x[0] - gam*x[1];
// 	}
// };

class FileObserver{
	std::ofstream &_fo;
public:
	FileObserver(std::ofstream& fo): _fo(fo) {}
	
	void operator()(State v, double t){
		_fo << t << '\t' << v[0] * cos(v[1]) << ' ' << v[0] * sin(v[1]) << ' ' << v[2] << ' ' << v[3] * cos(v[1]) - v[4] * sin(v[1]) << ' ' << v[3] * sin(v[1]) + v[3] * cos(v[1]) << ' ' << v[5] << '\n';
	}
};

class SpringForce{
	double _k;
public:
	/**
	 * Create isotropic elastic force attachment
	 * @param k elastic constant
	 */
	SpringForce(double k): _k(k) {}

	Vector3 operator()(State x, double /* t */ ){
		Vector3 f;
		f[0] = (-1 * _k) * x[0];
		return f;
	}
} force(1);

struct MagneticField{
	double B0;
	MagneticField(double B): B0(B) {}
	Vector3 operator()(Vector3 /* r */, double /* t */ ){
		Vector3 f;
		f[2] = B0;
		return f;
	}
} B(1);

struct ElectricField{
	double E0;
	ElectricField(double E): E0(E) {}
	Vector3 operator()(Vector3 r, double /* t */ ){
		Vector3 f;
		f[0] = E0 * cos(r[1]);
		f[1] = - E0 * sin(r[1]);
		return f;
	}
} E(10);

int main(){
	typedef Lorentz<NullForce, MagneticField, NullVectorField> System;
	System sys(1, 1, B, null_vector_field, null_force);

	EulerStepper<System, State, double> euler;
	RK46NL<System, State, double> rk46nl;
	Boris<System> boris;

	State x1 = {1, 0, 0, 0, 2, 1};
	State x2 = {1, 0, 0, 0, 2, 1};
	State x3 = {1, 0, 0, 0, 2, 1};

	std::ofstream euler_fo("xy_euler.dat");
	std::ofstream rk46nl_fo("xy_rk46.dat");
	std::ofstream boris_fo("xy_boris.dat");
	FileObserver obs_euler(euler_fo);
	FileObserver obs_rk46(rk46nl_fo);
	FileObserver obs_boris(boris_fo);

	size_t N = 100000;

	integrate(euler, sys, x1, 0.0, 0.01, N, obs_euler);

	auto start_rk46 = std::chrono::high_resolution_clock::now();
	integrate(rk46nl, sys, x2, 0.0, 0.01, N, obs_rk46, 9);
	auto end_rk46 = std::chrono::high_resolution_clock::now();

	auto start_boris = std::chrono::high_resolution_clock::now();
	integrate(boris, sys, x3, 0.0, 0.01, N, obs_boris, 9);
	auto end_boris = std::chrono::high_resolution_clock::now();

	std::cout << "Energy from Euler: " << mod(x1) << '\n';
	std::cout << "Runge-Kutta ran in " << std::chrono::duration_cast<std::chrono::milliseconds>(end_rk46 - start_rk46).count() << " ms\n";
	std::cout << "Energy : " << mod(x2) << '\n';
	std::cout << "Boris ran in " << std::chrono::duration_cast<std::chrono::milliseconds>(end_boris - start_boris).count() << " ms\n";
	std::cout << "Energy : " << mod(x3) << '\n';

	return 0;
}
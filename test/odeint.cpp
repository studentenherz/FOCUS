#include <iostream>
#include <fstream>

#include "odeint/integrator.hpp"
#include "odeint/stepper/euler.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "types/vector.hpp"

double gam = 0.0;

struct System
{
	System() {}
	void operator()(const Vector& x, Vector& dxdt, double t){
		dxdt[0] = x[1];
		dxdt[1] = -x[0] - gam*x[1];
	}
};

class FileObserver{
	std::ofstream &_fo;
public:
	FileObserver(std::ofstream& fo): _fo(fo) {}
	
	void operator()(Vector v, double t){
		_fo << t << '\t' << v << '\n';
	}
};


int main(){
	System sys;

	EulerStepper<System, Vector, double> euler;
	RK46NL<System, Vector, double> rk46nl;

	Vector x1(0, 1, 0);
	Vector x2 = x1;

	std::ofstream euler_fo("xy_euler.dat");
	std::ofstream rk46nl_fo("xy_rk46.dat");
	FileObserver obs_euler(euler_fo);
	FileObserver obs_rk46(rk46nl_fo);

	integrate(euler, sys, x1, 0.0, 0.01, 10000, obs_euler);
	integrate(rk46nl, sys, x2, 0.0, 0.01, 10000, obs_rk46);

	std::cout << "Energy from Euler: " << x1[0] * x1[0] + x1[1] * x1[1] << '\n';
	std::cout << "Energy from Runge-Kutta: " << x2[0] * x2[0] + x2[1] * x2[1] << '\n';

	return 0;
}
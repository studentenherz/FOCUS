#include <iostream>

#include "odeint/integrator.hpp"
#include "odeint/stepper/euler.hpp"
#include "types/vector.hpp"

double gam = 0.15;

struct System
{
	System() {}
	void operator()(const Vector& x, Vector& dxdt, double t){
		dxdt[0] = x[1];
		dxdt[1] = -x[0];
	}
};

void observer(Vector v, double t){
	std::cout << t << '\t' << v << '\n';
}


int main(){
	System sys;
	EulerStepper<System, Vector, Vector, double> euler;

	Vector x(0, 1, 0);

	integrate(euler, sys, x, 0.0, 0.01, 10000, observer);

	return 0;
}
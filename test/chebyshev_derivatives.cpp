#include <iostream>

#include "chebyshev.hpp"

int main(){
	int n = 10;
	double x = 0.67;
	double result[] = {0, 1, 2.68, 2.3868, -1.09558, -5.8131, -7.70409, -3.90572, 4.29079, 11.49, 11.7439};
	double e = 1e-6; // error tolerance

	Array<double> dT(n + 1);
	derivative_Chebyshev_T(n, x, dT);

	for (int i = 0; i <= n; i++)
		if(abs(dT[i] - result[i]) > e){
			std::cerr << "Got " << dT[i] << " expected " << result[i] << '\n';
			return 1;
		}

	return 0;
}

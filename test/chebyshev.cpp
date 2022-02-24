#include <iostream>
#include <iomanip>

#include "Chebyshev.hpp"

int main(int argc, char const *argv[]){
	double T[10], x = 0.452, n = 10;
	Chebyshev_T(n, x, T);
	
	std::cout << "Recurrence\tAnalytical\tDifference\n";
	std::cout << std::fixed;
	for (int i=0; i<=10; i++){
		double T_i_x = Chebyshev_T(i, x);
		std::cout << T[i] << '\t' << T_i_x << '\t' <<  T[i] - T_i_x << '\n';
	}

	return 0;
}

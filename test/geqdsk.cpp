#include "geqdsk.hpp"

int main(int argc, char* argv[]){
	if (argc < 2) return 0;

	Equilibrium eq = read_eqdsk(argv[1]);
	std::cout << eq.idnum << ' ' << eq.nx << ' ' << eq.ny << '\n';

	std:: cout << eq.bcentr << ' ' << eq.fpol[12] << ' ' << eq.psi(45, 67) << '\n';

	return 0;
}
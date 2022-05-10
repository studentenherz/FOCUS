#include "geqdsk.hpp"

int main(int argc, char* argv[]){
	if (argc < 2) return 0;

	Equilibrium eq = read_eqdsk(argv[1]);
	std::cout << eq.idnum << '\n';
	std::cout << eq.nx << '\n';
	std::cout << eq.ny << '\n';
	std::cout << eq.rdim << '\n';
	std::cout << eq.zdim << '\n';
	std::cout << eq.rcentr << '\n';
	std::cout << eq.rleft << '\n';
	std::cout << eq.zmid << '\n';
	std::cout << eq.rmagx << '\n';
	std::cout << eq.zmagx << '\n';
	std::cout << eq.simagx << '\n';
	std::cout << eq.sibdry << '\n';
	std::cout << eq.bcentr << '\n';
	std::cout << eq.cpasma << '\n';
	std::cout << eq.simagx << '\n';
	std::cout << eq.rmagx << '\n';
	std::cout << eq.zmagx << '\n';
	std::cout << eq.sibdry << '\n';

	return 0;
}
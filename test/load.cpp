#include <iostream>

#include "files.hpp"
#include "types.hpp"

int main(int argc, char const *argv[]){
	if (argc < 2) return 0;
	Matrix2D psi;
	if (load(argv[1], psi))
		std::cout << psi.size() << " x " << psi[0].size() << " matrix read\n";
	return 0;
}

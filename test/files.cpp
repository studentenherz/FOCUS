#include <iostream>

#include "files.hpp"
#include "types/matrix2D.hpp"

int main(int argc, char const *argv[]){
	if (argc < 2) return 0;
	Matrix2D<double> psi;
	if (!load(argv[1], psi))
		return 1;
	if (argc >= 3)
		dump(argv[2], psi);
	return 0;
}

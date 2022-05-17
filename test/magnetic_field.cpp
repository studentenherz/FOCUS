#include "files.hpp"
#include "magnetic_field.hpp"
#include "geqdsk.hpp"
#include "types/equilibrium.hpp"

int main(int argc, char* argv[]){
	if (argc < 2)
		return -1;
	
	Equilibrium eq = read_eqdsk(argv[1]);

	MagneticFieldMatrix B(eq, 26, 600);

	dump("Br.dat", B.Br, false);
	dump("Bz.dat", B.Bz, false);
	dump("Bt.dat", B.Bt, false);

	return 0;
}
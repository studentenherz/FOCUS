#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

#include "formats/matrix.hpp"
#include "formats/geqdsk.hpp"
#include "magnetic_field.hpp"
#include "types/equilibrium.hpp"
#include "types/vector.hpp"
#include "cxxopts.hpp"

bool make_dir(std::string dirname){
	std::string command = "mkdir -p ";
	command += dirname;
	return  (system(command.c_str()) == 0);
}

int main(int argc, char* argv[]){
	cxxopts::options options("Read from G-EQDSK file and return interpolated magnetif fields and poloidal flux");

	options.add_options()
		("g,geqdsk", "G-EQDSK file", cxxopts::value<std::string>())
		("o,output", "Output directory", cxxopts::value<std::string>())
		("h,help", "Display this help message")
	;

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			std::cout << options.help() << '\n';
			exit(0);
		}

		std::string ifilename = result["geqdsk"].as<std::string>();
		std::string odirname = result["output"].as<std::string>();

		if (!make_dir(odirname)){
			std::cerr << "Error creating directory " << odirname << "\n";
			exit(2);
		}

		Equilibrium eq = read_geqdsk(ifilename);
		FineEquilibrium fineq(eq, 26);
		MagneticFieldMatrix B_matrix(eq, 26, 600);

		std::ofstream interpo_lims(odirname + "/interp_lims.dat");
		interpo_lims << B_matrix.r_min * eq.rdim << ' ' << B_matrix.r_max * eq.rdim << ' ' << B_matrix.z_min * eq.rdim << ' ' << B_matrix.z_max * eq.rdim << '\n';

		interpo_lims.close();

		dump(odirname + "/Br.dat", B_matrix.Br, false);
		dump(odirname + "/Bt.dat", B_matrix.Bt, false);
		dump(odirname + "/Bz.dat", B_matrix.Bz, false);
		dump(odirname + "/Psi.dat", B_matrix.psi, false);


		std::ofstream fobdry(odirname + "/bdry.dat");
		for (size_t i = 0; i < eq.nbdry; i++)
			fobdry << eq.rbdry[i] << ' ' << eq.zbdry[i] << '\n';
		fobdry.close();

		std::ofstream folim(odirname + "/lim.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			folim << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';
		folim.close();
	}
	catch (cxxopts::option_error const& e){
		std::cerr << e.what() << '\n';
		exit(1);
	}
	

	return 0;
}
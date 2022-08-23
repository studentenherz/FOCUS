#include <fstream>
#include <string>

#include "geqdsk.hpp"
#include "files.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){
	
	cxxopts::options options("geqdsk", "Test geqdsk read");

	options.add_options()
		("file", "", cxxopts::value<std::string>())
		("h,help", "Show this help message");

	options.positional_help("<G-EQDSK input file>");
	options.parse_positional({"file"});

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		Equilibrium eq = read_geqdsk(result["file"].as<std::string>().c_str());
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

		for (size_t i = 0; i < eq.nx; i++)
			std::cout << eq.fpol[i] << '\t' << eq.pres[i] << '\t' << eq.qpsi[i] << '\n';

		dump("psi_from_geqdsk.dat", eq.psi, false);
		
		std::ofstream fo1("psi.dat");
		for (size_t j = 0; j < eq.ny; j++)
			for (size_t i = 0; i < eq.nx; i++){
				double r = eq.rleft + i * eq.rdim / eq.nx;
				double z = eq.zmid - eq.zdim / 2 + j * eq.zdim / eq.ny;
				fo1 << r << ' ' << z << ' ' << eq.psi(i, j) << '\n';
			}


		// Output limit
		std::ofstream fo("limit.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			fo << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';

		std::ofstream fo2("boundary.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			fo2 << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';
		

		return 0;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
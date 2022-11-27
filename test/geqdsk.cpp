#include <fstream>
#include <string>

#include "formats/geqdsk.hpp"
#include "formats/matrix.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){
	
	cxxopts::options options("geqdsk", "Get poloidal flux and other information out of G-EQDSK files");

	options.add_options()
		("g,geqdsk", "G-EQDSK input file", cxxopts::value<std::string>())
		("o,output", "Output directory", cxxopts::value<std::string>())
		("h,help", "Show this help message");

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		Equilibrium eq = read_geqdsk(result["geqdsk"].as<std::string>());

		std::cout << "idnum " << eq.idnum   << '\n';
		std::cout << "nx    " << eq.nx      << '\n';
		std::cout << "ny    " << eq.ny      << '\n';
		std::cout << "rdim  " << eq.rdim    << '\n';
		std::cout << "zdim  " << eq.zdim    << '\n';
		std::cout << "rcentr" << eq.rcentr  << '\n';
		std::cout << "rleft " << eq.rleft   << '\n';
		std::cout << "zmid  " << eq.zmid    << '\n';
		std::cout << "rmagx " << eq.rmagx   << '\n';
		std::cout << "zmagx " << eq.zmagx   << '\n';
		std::cout << "simagx " << eq.simagx  << '\n';
		std::cout << "sibdry " << eq.sibdry  << '\n';
		std::cout << "bcentr " << eq.bcentr  << '\n';
		std::cout << "cpasma " << eq.cpasma  << '\n';
		std::cout << "simagx " << eq.simagx  << '\n';
		std::cout << "rmagx " << eq.rmagx   << '\n';
		std::cout << "zmagx " << eq.zmagx   << '\n';
		std::cout << "sibdry " << eq.sibdry  << '\n';

		std::string odirname = result["output"].as<std::string>();

		std::ofstream psilims(odirname + "/psilims.dat");
		psilims << eq.rleft << ' ' << eq.rleft + eq.rdim << ' ' << eq.zmid - 0.5 * eq.zdim << ' ' << eq.zmid + 0.5 * eq.zdim;
		psilims.close();

		dump(odirname + "/psi_from_geqdsk.dat", eq.psi, false);

		// Output limit
		std::ofstream fo(odirname + "/lim.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			fo << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';

		std::ofstream fo2(odirname + "/bdry.dat");
		for (size_t i = 0; i < eq.nbdry; i++)
			fo2 << eq.rbdry[i] << ' ' << eq.zbdry[i] << '\n';

		return 0;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
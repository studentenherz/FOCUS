#include <iostream>
#include <sstream>

#include "collisions.hpp"
#include "cxxopts.hpp"
#include "formats/geqdsk.hpp"
#include "formats/input_gacode.hpp"
#include "formats/particle_states.hpp"
#include "lorentz.hpp"
#include "magnetic_field.hpp"
#include "odeint/integrator.hpp"
#include "odeint/stepper/collision_stepper.hpp"
#include "odeint/stepper/rk46_nl.hpp"
#include "types/array.hpp"
#include "types/equilibrium.hpp"
#include "types/particle.hpp"
#include "types/plasma.hpp"
#include "types/vector.hpp"

int main(int argc, char* argv[]){
	std::ostringstream oss;
	oss << "\n\tFOCUS: Full Orbit CUda Solver (" << GIT_VERSION << ")\n";
	std::string hello_message = oss.str();

	cxxopts::options options(argv[0], hello_message);

	options
		.set_width(90)
    .set_tab_expansion()
		.add_options()
		("g,geqdsk", "G-EQDSK input file", cxxopts::value<std::string>(), "ifile")
		("i,input-gacode", "input.gacode input file", cxxopts::value<std::string>(), "ifile")
		("x0,initial-states", "initial states input file", cxxopts::value<std::string>(), "ifile")
		("o,output", "Output file for particles states", cxxopts::value<std::string>(), "ofile")
		("b", "Dump magnetic fields and poloidal flux")
		("t0", "Initial time", cxxopts::value<double>()->default_value("0.0"))
		("dt", "Time steps", cxxopts::value<double>()->default_value("0.0001"))
		("N", "Number of steps of integration", cxxopts::value<size_t>()->default_value("1000000"))
		("h,help", "Show this help message");

	try{
		auto result = options.parse(argc, argv);

		if(result.count("help")){
			std::cout << options.help() << std::endl;
			return EXIT_SUCCESS;
		}

		std::string geqdsk = result["geqdsk"].as<std::string>();
		std::string input_gacode = result["input-gacode"].as<std::string>();

		// Read equilibrium
		std::cout << "Reading G-EQDSK input file " << std::flush;
		Equilibrium eq = read_geqdsk(geqdsk);
		std::cout << "...done\n";

		// Plasma particles
		std::cout << "Reading input.gacode input file " << std::flush;
		Plasma plasma = read_input_gacode(input_gacode);
		std::cout << "...done\n";

		// Refined grid and B
		std::cout << "Computing Chebyshev expansion " << std::flush;
		MagneticFieldMatrix B_matrix(eq, 26, 600);
		std::cout << "...done\n";

		// 

		return EXIT_SUCCESS;
	}catch(cxxopts::option_error const& e){
		std::cerr << '\n' << e.what()
		<< "\n\nRun ‘" << argv[0] << " --help‘ in order to see the required options.\n";
		return EXIT_FAILURE;
	}
}
#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "formats/particle_states.hpp"
#include "types/array.hpp"
#include "types/vector.hpp"
#include "types/particle.hpp"

int main(int argc, char* argv[]){
	cxxopts::options options("input_states", "Test reading input states");

	options.add_options()
		("input", "", cxxopts::value<std::string>())
		("output", "", cxxopts::value<std::string>()->default_value("states.out"));

	options.positional_help("<states input file> [<states output file>]");
	options.parse_positional("input", "output");

	try{
		auto result = options.parse(argc, argv);

		Array<State> states = load_states(result["input"].as<std::string>());

		Particle part(1, 1.2);

		dump_states(result["output"].as<std::string>(), states, part);

		return EXIT_SUCCESS;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}
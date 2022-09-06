#include <iostream>
#include <string>

#include "cxxopts.hpp"
#include "formats/input_states.hpp"
#include "types/array.hpp"
#include "types/vector.hpp"

int main(int argc, char* argv[]){
	cxxopts::options options("input_states", "Test reading input states");

	options.add_options()
		("input", "", cxxopts::value<std::string>());

	options.positional_help("<states input file>");
	options.parse_positional("input");

	try{
		auto result = options.parse(argc, argv);

		Array<State> states = load_states(result["input"].as<std::string>());

		for(size_t i =0 ; i < states.size(); i++)
			std::cout << states[i] << '\n';

		return EXIT_SUCCESS;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}
#include <string>

#include "formats/input_gacode.hpp"
#include "cxxopts.hpp"
#include "types/plasma.hpp"

int main(int argc, char* argv[]){
	cxxopts::options options("input_gacode", "Test reading input.gacode files");

	options.add_options()
		("file", "", cxxopts::value<std::string>());

	options.parse_positional({"file"});

	try{
		auto result = options.parse(argc, argv);
		
		std::vector<std::string> species_identifiers;
		Plasma plasma = read_input_gacode(result["file"].as<std::string>(), species_identifiers);

		std::cout << "nexp " << plasma.nexp << '\n';
		std::cout << "nion " << plasma.nion << '\n';
		std::cout << "shot " << plasma.shot << '\n'; 
		std::cout << "masse " << plasma.masse << '\n'; 
		std::cout << "ze " << plasma.ze << '\n';
		 
		std::cout << "species identifiers\n";
		for(size_t i = 0; i< plasma.nion; i++)
			std::cout << species_identifiers[i] << ' ';
		 
		std::cout << "\nmass\n";
		for(size_t i = 0; i< plasma.nion; i++)
			std::cout << plasma.mass[i] << '\n';

		std::cout << "z\n";
		for(size_t i = 0; i< plasma.nion; i++)
			std::cout << plasma.z[i] << '\n';

		std::cout << "psi\n";
		for(size_t i = 0; i< plasma.nexp; i++)
			std::cout << plasma.polflux[i] << '\n';

		std::cout << "ne\n";
		for(size_t i = 0; i< plasma.nexp; i++)
			std::cout << plasma.ne[i] << '\n';

		std::cout << "te\n";
		for(size_t i = 0; i< plasma.nexp; i++)
			std::cout << plasma.te[i] << '\n';

		std::cout << "ni\n";
		for (size_t i = 0; i < plasma.nexp; i++){
				for (size_t ion = 0; ion < plasma.nion; ion++){
					std::cout << plasma.ni(ion, i) << '\t';
				}
			std::cout << '\n';
		}

		std::cout << "ti\n";
		for (size_t i = 0; i < plasma.nexp; i++){
				for (size_t ion = 0; ion < plasma.nion; ion++){
					std::cout << plasma.ti(ion, i) << '\t';
				}
			std::cout << '\n';
		}

		return EXIT_SUCCESS;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << '\n';
		return EXIT_FAILURE;
	}
}
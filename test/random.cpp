#include <iostream>
#include <fstream>
#include <string>

#include "random.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){
	cxxopts::Options options("random", "Generate random distributions from seed using ran2");

	typedef unsigned long long ull;

	options.add_options()
	("s,seed", "Seed of random", cxxopts::value<ull>())
	("h,help", "Show this help message")
	("uniform_ofile", "Uniform distribution output file", cxxopts::value<std::string>()->default_value("uniform.dat"));
	("normal_ofile", "Normal distribution output file", cxxopts::value<std::string>()->default_value("normal.dat"));

	try{

		auto result = options.parse(argc, argv);

		if(result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		ull seed = result["seed"].as<ull>();
		
		Ran2 ran(seed);
		std::ofstream uniform(result["uniform_ofile"].as<std::string>());
		for(int i =0 ; i<10000; i++)
			uniform << ran.doub() << '\n';

		NormalRand norm(seed, 1.0, 5.0);
		std::ofstream normal(result["normal_ofile"].as<std::string>());
		for(int i =0 ; i<10000; i++)
			normal << norm() << '\n';

		return 0;
	} catch(cxxopts::OptionException const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
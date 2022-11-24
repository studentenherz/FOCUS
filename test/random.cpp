#include <iostream>
#include <fstream>
#include <string>

#include "random.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){
	cxxopts::options options("random", "Generate random distributions from seed using ran2");

	typedef unsigned long long ull;

	options.add_options()
	("s,seed", "Seed of random", cxxopts::value<ull>(), "n")
	("u,uniform_ofile", "Uniform distribution output file", cxxopts::value<std::string>()->default_value("uniform.dat"), "file")
	("n,normal_ofile", "Normal distribution output file", cxxopts::value<std::string>()->default_value("normal.dat"), "file")
	("h,help", "Show this help message");

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
			uniform << ran.uniform() << '\n';

		NormalRand norm(seed, 1.0, 5.0);
		std::ofstream normal(result["normal_ofile"].as<std::string>());
		for(int i =0 ; i<10000; i++)
			normal << norm() << '\n';

		return 0;
	} catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
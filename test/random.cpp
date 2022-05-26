#include <iostream>
#include <fstream>

#include "random.hpp"

int main(int argc, char* argv[]){
	if(argc < 2){
		std::cout << "Usage:\nrandom <seed>\n";
		return 0;
	}

	unsigned long long seed = std::atoll(argv[1]);
	
	Ran2 ran(seed);
	std::ofstream uniform("uniform.dat");
	for(int i =0 ; i<10000; i++)
		uniform << ran.doub() << '\n';

	NormalRand norm(seed, 1.0, 5.0);
	std::ofstream normal("normal.dat");
	for(int i =0 ; i<10000; i++)
		normal << norm() << '\n';


	return 0;
}
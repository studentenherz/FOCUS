#include <iostream>

#include "types/vector.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){
	using namespace std;

	cxxopts::options options("vector", "A simple test of vectors in this library");

	options.add_options()
		("b,bfirst", "First coordinate of b", cxxopts::value<double>()->default_value("0"))
		("h,help", "Show this help");

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			cout << options.help() << endl;
			exit(0);
		}

		
		Vector3 zero;
		Vector3 a = {1.2, 4.0, -9}, b;
		b = {result["bfirst"].as<double>(), 1, 3};

		cout << "zero = " << zero << '\n';
		cout << "a = " << a << '\n';
		cout << "b = " << b << '\n';

		cout << "2 * a = " << 2 * a << '\n';
		cout << "a * (-2) = " << a * (-2) << '\n';
		cout << "a / 2 = " << a / 2 << '\n';

		cout << "a + b = " << a + b << '\n';
		cout << "a - b = " << a - b << '\n';
		cout << "dot(a, b) = " << dot(a, b) << '\n';
		cout << "cross(a, b) = " << cross(a, b) << '\n';
		cout << "mod(a) = " << mod(a) << '\n';
		return 0;
	}
	catch(cxxopts::option_error const& e){
		cerr << e.what() << endl;
		return 1;
	}

}
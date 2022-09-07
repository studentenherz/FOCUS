#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "types/array.hpp"
#include "types/vector.hpp"
#include "magnetic_field.hpp"
#include "types/particle.hpp"
#include "formats/regex_tokenizer.hpp"
#include "formats/geqdsk.hpp"
#include "formats/particle_states.hpp"
#include "util.hpp"
#include "random.hpp"
#include "cxxopts.hpp"



//   R(cm)         Z(cm)         vpll/v        E(eV)         zeta(deg)
State pitch_E_2_state(double R, double z, double vpll_v, double E, double theta, Particle part, MagneticFieldFromMatrix& B, double a, Ran2& ran){
	double theta_rad = theta * pi / 180.0;
	double r_m = R / 100.0;
	double z_m = z / 100.0;

	double v_mod = 4.39284e5 * sqrt(E/part.m);
  Vector3 r = {r_m / a, theta_rad, z_m / a}; // Normalized
	Vector3 b = B(r, 0.0);

  Vector3 e_ll = b / mod(b);
  // First perpendicular (T) unit vector
  Vector3 e_x = {1, 0, 0};
  Vector3 e_T_1 = cross(e_ll, e_x);
  if (mod(e_T_1) == 0){
    Vector3 e_y = {0, 1, 0};
    e_T_1 = cross(e_ll, e_y);
  }
  e_T_1 = e_T_1 / mod(e_T_1);

  // Second perpendicular (T) unit vector
	Vector3 e_T_2 = cross(e_ll, e_T_1);

  double v_ll = vpll_v * v_mod;
  double v_T = v_mod * sqrt(1 - sqr(vpll_v));
  
  double phi = ran.random(0, two_pi);

  double v_T_1 = v_T * cos(phi);
  double v_T_2 = v_T * sin(phi);

  Vector3 v = v_ll * e_ll + v_T_1 * e_T_1 + v_T_2 * e_T_2;

	State state{r_m, theta_rad, z_m, v[0], v[1], v[2]};

	if (hasnan(state)){
		std::cout << state << '\n';
		std::cout << "v_mod = " << v_mod << '\n';
		std::cout << "r = " << r << '\n';
		std::cout << "b = " << b << '\n';
		std::cout << "e_ll = " << e_ll << '\n';
		std::cout << "e_T_1 = " << e_T_1 << '\n';
		std::cout << "e_T_2 = " << e_T_2 << '\n';
		std::cout << "v_ll = " << v_ll << '\n';
		std::cout << "v_T = " << v_T << '\n';
		std::cout << "phi = " << phi << '\n';
		std::cout << "v_T_1 = " << v_T_1 << '\n';
		std::cout << "v_T_2 = " << v_T_2 << '\n';
		std::cout << "v = " << v << '\n';
		exit(EXIT_FAILURE);
	}

  return state;
}

Array<State> load_states_from_pitch_E(std::string filename, Particle part, MagneticFieldFromMatrix& B, double a, Ran2& ran){
	std::ifstream fi(filename);
	if(!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return Array<State>(0);
	}

	std::vector<State> v_states;

	double R, z, vpll_v, E, theta;
	std::string line;

	while	(std::getline(fi, line) && line.find("<start-of-data>") == line.npos) {}

	while	(std::getline(fi, line)){
		std::istringstream iss(line);
		iss >> R >> z >> vpll_v >> E >> theta;
    State stte = pitch_E_2_state(R, z, vpll_v, E, theta, part, B, a, ran);
		v_states.push_back(stte);
	}

	Array<State> states(v_states.size());
	for (size_t i = 0; i < states.size(); i++)
		states[i] = v_states[i];

	return states;
}


int main(int argc, char* argv[]){
  cxxopts::options options(argv[0], "\n\tConvert transp pith-E input file into inital states input file\n");

  options.add_options()
    ("g,geqdsk", "G-EQDSK input file", cxxopts::value<std::string>())
    ("p,pitch-e", "pitch-E input file", cxxopts::value<std::string>())
    ("o,output", "Output stetes file", cxxopts::value<std::string>())
    ("m", "Particle mass [Da]", cxxopts::value<double>()->default_value("2.0147294"))
    ("Z", "Particle Z", cxxopts::value<double>()->default_value("1"))
    ("h,help", "Display this help message");

  try{
    auto result = options.parse(argc, argv);

    if (result.count("help")){
      std::cout << options.help() << std::endl;
      return EXIT_SUCCESS;
    }

    std::string geqdsk = result["geqdsk"].as<std::string>();
    std::string pitch_e = result["pitch-e"].as<std::string>();
    std::string output = result["output"].as<std::string>();

    Equilibrium eq = read_geqdsk(geqdsk);
    MagneticFieldMatrix B_matrix(eq, 26, 600);
    MagneticFieldFromMatrix B(B_matrix, eq.bcentr);

    Particle part(result["Z"].as<double>(), result["m"].as<double>());
    Ran2 ran(12);

    Array<State> states = load_states_from_pitch_E(pitch_e, part, B, eq.rdim, ran);

    dump_states(output, states, part);

    return EXIT_SUCCESS;
  }catch(cxxopts::option_error const& e){
    std::cerr << e.what() << "\n\nRun ‘" << argv[0] << " --help‘ in order to see the required options.\n";
    return EXIT_FAILURE;
  }

}
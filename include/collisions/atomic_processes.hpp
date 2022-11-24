#if !defined(FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP)
#define FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP

#include <string>
#include <fstream>
#include <sstream>

#include "types/matrix_2d.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/plasma.hpp"
#include "formats/matrix.hpp"
#include "types/vector.hpp"
#include "interpolations.hpp"

enum AtomicProcessType {
	ChargeExchange, 
	DeExcitation,		// The probability of this one is calculated in a different fashion from the others
	Excitation,
	ElectronCapture,
	Ionization,
	Other
};

class AtomicProcess{
	/** 
	 * Data for the process, it might contain a single scalar in case of de-excitation, 
	 * an array or a 	matrix in case of 1D or 2D interpolation is needed respectively.
	 */
	Matrix2D<double> M;			
	
	Array<double> E;				// Energies for interpolation
	Array<double> T;				// Temperatures for interpolation

	size_t other_index; 		// Index of the other reactant
	
	AtomicProcessType type;	// Atomic process type

	int q;									// Charge for which the process applies
	ulong n;								// Quantum number for which the process applies

	int next_q;						// Resulting charge if the process does occur
	ulong next_n;					// Resulting quantum number if the process occurs
public:

	AtomicProcess() {}

	/**
	 * Constructor for AtomicProcess 
	 * (not a C++ class constructor but a method that populates the object)
	 * 
	 * @param input_filename file with the process rates information
	 * @param other_idex index in the Plasma of the other reactant species
	 * @param q charge for which the process applies
	 * @param n quantum number for which the process applies
	 * @param next_q resutling charge after the process
	 * @param next_n resulting state after the process
	 * @param type type of atomic process
	 */
	#ifdef __CUDACC__
	__host__
	#endif
	void constructor(std::string input_filename, size_t _other_index, int _q, ulong _n, int _next_q, ulong _next_n, AtomicProcessType _type){
		/**
		 * The input file for each atomic process should be structured in the following way:
		 * 
		 * Nx Ny
		 * [Array E] (dim: Ny)
		 * [Array T] (dim: Nx)
		 * Matrix M (dim: Nx x Ny)
		 * 
		 * Nx, Ny are the dimension of the data
		 * 
		 * If the data is only a sacalar, Nx = Ny = 1 and the E, T are not given
		 * If Nx > 1 but Ny = 1, E is given, and it means a dependency with E to be interpolated
		 * If Nx > 1 and Ny > 1 E and T are given for bilinear interpolation 
		*/

		other_index = _other_index;
		q = _q;
		n = _n;
		next_q = _next_q;
		next_n = _next_n;
		type = _type;

		std::ifstream fi(input_filename);
		if (!fi.is_open()){
			std::cerr << "Unable to open file " << input_filename << "\n";
			exit(1);
		}

		// Read dimensions of the matrix
		size_t Nx, Ny;
		fi >> Nx >> Ny;

		if (Ny > 1){ // Read E
			E.resize(Ny);
			for (size_t i = 0; i < Ny; i++)
				fi >> E[i];
		}

		if (Nx > 1){ // Read T
			T.resize(Nx);
			for (size_t i = 0; i < Nx; i++)
				fi >> T[i];
		}

		// Read M
		M.reshape(Nx, Ny);
		for (size_t j = 0; j < Ny; j++)
			for (size_t i = 0; i < Nx; i++)
				fi >> M(i, j);

	}

	#ifdef __CUDACC__
	__host__
	cudaError_t construct_in_host_for_device(AtomicProcess& other){
		other_index = other.other_index;
		type = other.type;
		q = other.q;
		n = other.n;
		next_q = other.next_q;
		next_n = other.next_n;

		propagateCUDAErr( M.construct_in_host_for_device(other.M) );
		propagateCUDAErr( E.construct_in_host_for_device(other.E) );
		propagateCUDAErr( T.construct_in_host_for_device(other.T) );

		return cudaSuccess;
	}
	#endif

	/**
	 * Calculate the probability for this particular atomic process
	 * 
	 * @param part test particle
	 * @param energy energy of the test particle
	 * @param psi poloidal flux at the position of the particle
	 * @param plasma Plasma through which the particle is moving
	 * @param dt process calculation delta time
	 * @return Probability of the process (<< 1)
	*/
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double P(const Particle& part, double energy, double psi, Plasma& plasma, double t, double dt, double n_conversion_factor_4_APs = 1){
		// This process only applies to particles in state q, n
		if (part.q != q || part.n != n) return 0;

		if (type == AtomicProcessType::DeExcitation)
			return 1.0 - exp(- part.t / M(0, 0));

		// Density of plasma species
		double n = plasma[other_index].n(psi, t) * n_conversion_factor_4_APs;

		// Linear interpolation
		if (T.size() == 0)
			return n * lagrange_interpolation_3(energy, E, M, 0) * dt;
		
		// Bilinear interpolation
		double temp = plasma[other_index].T(psi, t);
		return n * four_point_formula(temp, energy, T, E, M) * dt;
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void apply(Particle& part){
		part.q = next_q;
		part.n = next_n;
		part.t = 0;
	}
};

/**
 * @brief Loads the atomic processes for a specific species
 * 
 * Inside the directory there should be a file names .index for ease of 
 * reading, .index is structured like this:
 * 
 * nfiles 
 * filename_1
 * # Comments start with #
 * filename_2
 * ...
 * filename_n
 * 
 * nfiles is the number filenames to come
 * 
 * Each filename afterwards represents a file inside the same directory with 
 * the needed information about the atomic process.
 * 
 * Each filename is in the following format
 * 
 * [type of process (Io|CE|Ex|DE)]_[identifier of particle type]_[identifier of the other reactant particle type]_[q]_[n]_[next_q]_[next_n]
 * 
 * For example:
 * 
 * 		Ex_D_e_0_2_0_3
 * 
 * would have the data for the process of excitation of a neutral Deuterium by electrons from n = 2 to n = 3 (both states q = 0)
 * 
 * The type of process (for now) is only relevant when the process is a DeExcitation, as it is calculated
 * in a different fashion than the other processes.
 * 
 * @param plasma Plasma object
 * @param directory Directory from which to retrieve the information
 * 
 * @returns Array of pointers to atomic processes
 */
Array<AtomicProcess> load_atomic_processes(std::vector<std::string> species_identifiers, std::string directory){
	std::string index_filename = directory + "/.index";
	std::ifstream index_file (index_filename);
	if (!index_file.is_open()){
		std::cerr << "Error loading atomic processes, couldn't open " << index_filename << '\n';
		exit(1);
	}
	size_t nfiles;
	index_file >> nfiles;

	std::cout << "Loading " << nfiles << " atomic processes\nPlasma species: ";
	for (auto id: species_identifiers)
		std::cout << id << ' ';
	std::cout << '\n'; 

	Array<AtomicProcess> processes(nfiles);

	std::string filename;
	AtomicProcessType type;
	size_t other_index = 0;
	ulong q, n, next_q, next_n;

	
	size_t count = 0;
	for (size_t i = 0; i < nfiles; i++){
		std::getline(index_file, filename);
		if (filename.empty() || filename[0] == '#'){ // Ignore comments and empty lines
			i--;
			continue;
		}

		std::stringstream ss(filename);
		std::string token;
		
		// Type
		std::getline(ss, token, '_');
		if (token == "DE")
			type = AtomicProcessType::DeExcitation;
		else
			type = AtomicProcessType::Other;

		// Test particle species. Won't use for now
		std::getline(ss, token, '_');

		// Other particle species
		std::getline(ss, token, '_');
		if (type != AtomicProcessType::DeExcitation){ // This does not matter for de-excitation
			if (token == "e") other_index = species_identifiers.size();
			else{
				other_index = 0;
				while (other_index < species_identifiers.size() && species_identifiers[other_index] != token) other_index++;
				if (other_index == species_identifiers.size()) // Not present in this plasma, next file
					continue;
			}
		}

		std::getline(ss, token, '_');
		q = std::stoi(token);

		std::getline(ss, token, '_');
		n = std::stoul(token);
		
		std::getline(ss, token, '_');
		next_q = std::stoi(token);
		
		std::getline(ss, token, '_');
		next_n = std::stoul(token);

		processes[i].constructor(directory + "/" + filename, other_index, q, n, next_q, next_n, type);

		std::cout << ++count << " Loaded " << filename << '\n';
	}

	if (count != nfiles){
		std::cerr << nfiles << " atomic processes files were spected but only got " << count << '\n';
		exit(1);
	}

	return processes;
}

template <typename RandomGenerator_t, typename MagneticField_t>
class AtomicProcessesHandler{
	Array<AtomicProcess>& processes;
	RandomGenerator_t& ran_gen;
	Plasma& plasma;
	double energy_conversion_factor;
	double n_conversion_factor_4_APs;
	double Omega;
public:
	MagneticField_t& B;
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	AtomicProcessesHandler(Array<AtomicProcess>& processes, RandomGenerator_t& ran_gen, MagneticField_t& B, Plasma& plasma, double Omega, double energy_conversion_factor, double n_conversion_factor_4_APs): processes(processes), ran_gen(ran_gen), B(B), plasma(plasma), Omega(Omega), energy_conversion_factor(energy_conversion_factor), n_conversion_factor_4_APs(n_conversion_factor_4_APs) {}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	void operator()(Particle& part, State& x, double t, double dt){
		dt /= Omega; // Time in seconds
		part.t += dt;

		Vector3 r = get_position(x);
		Vector3 v = get_velocity(x);

		double energy = part.m * dot(v, v) * energy_conversion_factor;
		double psi = B.psi(r, t);

		Array<double> p(processes.size());
		for (size_t i = 0; i < processes.size(); i++)
			p[i] = processes[i].P(part, energy, psi, plasma, t, dt, n_conversion_factor_4_APs);
		
		for (size_t i = 1; i < processes.size(); i++)
			p[i] += p[i - 1];

		double ran = ran_gen.uniform();

		for (size_t i = 0 ; i < processes.size(); i++)
			if (ran < p[i]){
				processes[i].apply(part);
				return;
			}
	}
};

#endif // FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP

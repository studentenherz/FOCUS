#if !defined(FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP)
#define FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP

#include <string>
#include <fstream>

#include "types/matrix_2d.hpp"
#include "types/array.hpp"
#include "types/particle.hpp"
#include "types/plasma.hpp"
#include "formats/matrix.hpp"
#include "types/vector.hpp"
#include "interpolations.hpp"

enum AtomicProcessType {
	Ionization,
	ChargeExchange, 
	DeExcitation		// The probability of this one is calculated in a different fashion from the others
};

class AtomicProcess{
	Matrix2D<double> M;			// Data for the process, it might contain a single scalar in case of de-excitation,  
													// an array or a matrix in case of 1D or 2D interpolation is needed respectively.
	
	Array<double> E;				// Energies for interpolation
	Array<double> T;				// Temperatures for interpolation

	size_t other_index; 		// Index of the other reactant
	
	AtomicProcessType type;	// Atomic process type

	uint q;									// Charge for which the process applies
	uint n;									// Quantum number for which the process applies

	uint next_q;						// Resulting charge if the process does occur
	uint next_n;						// Resulting quantum number if the process occurs
public:

	#ifdef __CUDACC__
	__host__
	#endif
	AtomicProcess(std::string input_filename, size_t other_index, uint q, uint n, uint next_q, uint next_n): other_index(other_index), q(q), n(n), next_q(next_q), next_n(next_n) {
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
		std::ifstream fi(input_filename);

		// Read dimensions of the matrix
		size_t Nx, Ny;
		fi >> Nx >> Ny;

		if (Nx > 1){ // Read E
			E.reshape(Nx);
			for (size_t i = 0; i < Nx; i++)
				fi >> E[i];
		}

		if (Ny > 1){ // Read T
			T.reshape(Ny);
			for (size_t i = 0; i < Ny; i++)
				fi >> T[i];
		}

		// Read M
		M.reshape(Nx, Ny);
		for (size_t j = 0; j < Ny; j++)
			for (size_t i = 0; i < Nx; i++)
				fi >> M(i, j);
	}

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
	double P(const Particle& part, double energy, double psi, const Plasma& plasma, double dt){
		// This process only applies to particles in state q, n
		if (part.q != q || part.n != n) return 0;

		if (type == AtomicProcessType::DeExcitation)
			return 1 - exp(- part.t / M(0, 0));

		// Density of plasma species
		double n = plasma[other_index].n(psi);

		// Linear interpolation
		if (T.size() == 0)
			return n * lagrange_interpolation_3(energy, E, M, 0) * dt;
		
		// Bilinear interpolation
		double temp = plasma[other_index].T(psi);
		return n * four_point_formula(t, energy, T, E, M) * dt;
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

#endif // FOCUS_INCLUDE_COLLISIONS_ATOMIC_PROCESSES_HPP

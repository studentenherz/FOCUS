#if !defined(FOCUS_INCLUDE_MAGNETIC_FIELD_HPP)
#define FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

#include "chebyshev.hpp"
#include "types/equilibrium.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "types/vector.hpp"

struct MagneticFieldMatrix{
	Matrix2D<double> Br;	// B_r
	Matrix2D<double> Bt;	// B_t
	Matrix2D<double> Bz;	// B_z

	/**
	 * Calculate magnetic fields from n'th order
	 * Chebyshev expansion of psi
	 * @param eq Equilibrium from EFIT's G-EQDSK
	 * @param n order of Chebyshev expansion
	 * @param N size of the refined grid
	 * @param sign if true (default) the negative sign on the poloidal
	 * fields definition goes to the radial component
	 */
	MagneticFieldMatrix(Equilibrium &eq, size_t n, size_t N, bool sign = true): Br(N, N), Bt(N, N), Bz(N, N) {
		double r_min = eq.rleft;
		double r_max = eq.rleft + eq.rdim;
		double z_min = eq.zmid - 0.5 * eq.zdim;
		double z_max = eq.zmid + 0.5 * eq.zdim;

		ScalarField psi(&eq.psi, r_min, r_max, z_min, z_max);
		ChebyshevExpansion ch(n, psi, r_min, r_max, z_min, z_max);

		double D_psi = (eq.sibdry - eq.simagx);
		double d_psi = D_psi / eq.nx;

		for(size_t i = 0; i<N; i++){
			double r = r_min + eq.rdim * i / (N - 1);
			for(size_t j = 0; j<N; j++){
				double z = z_min + eq.zdim * j / (N - 1);

				// From \Psi definition
				Br(i, j) =  (sign ? -1 : 1) * ch.dy(r, z) / r;
				Bz(i, j) =  (sign ? 1 : -1) * ch.dx(r, z) / r;

				double psi_here = ch(r, z);
				size_t index = std::floor((psi_here - eq.simagx) / D_psi);
				// Keep index inside boundaries
				if (index == 0) index = 1;
				if (index >= eq.nx - 1) index = eq.nx - 1;

				// Lagrange interpolation with the 3 closest points
				// https://en.wikipedia.org/wiki/Lagrange_polynomial
				Vector3 x, y, l;

				// Get the points
				for(size_t k = 0; k < 3; k++){
					x[k] = eq.simagx + (index + k - 1) * d_psi;
					y[k] = eq.fpol[index + k  - 1];
				}

				// Basis polynomials evaluated at psi_here
				for(size_t k = 0; k < 3; k++)
					l[k] = (psi_here - x[(k + 2) % 3])/(x[k] - x[(k + 2) % 3]) * (psi_here - x[(k + 1) % 3])/(x[k] - x[(k + 1) % 3]);
				
				// Linear combination
				double F = dot(y, l);

				// From F definition
				Bt(i, j) = F / r;
			}
		}
	}
};

#endif // FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

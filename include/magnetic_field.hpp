#if !defined(FOCUS_INCLUDE_MAGNETIC_FIELD_HPP)
#define FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

#include "chebyshev.hpp"
#include "interpolations.hpp"
#include "types/equilibrium.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "types/vector.hpp"

struct MagneticFieldMatrix{
	Matrix2D<double> Br;	// B_r
	Matrix2D<double> Bt;	// B_t
	Matrix2D<double> Bz;	// B_z

	double r_min, r_max; // r limits
	double z_min, z_max; // z limits

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

		// Dimensionless limits of the matrixes 
		double mr_min = eq.rleft / eq.rdim;
		double mr_max = eq.rleft / eq.rdim + 1;
		double mz_min = (eq.zmid - 0.5 * eq.zdim) / eq.rdim;
		double mz_max = (eq.zmid + 0.5 * eq.zdim) / eq.rdim;

		ScalarField psi(eq.psi, mr_min, mr_max, mz_min, mz_max);
		
		// Dimensionless limits of plasma boundaries
		r_min = min(eq.rbdry) / eq.rdim;
		r_max = max(eq.rbdry) / eq.rdim;
		z_min = min(eq.zbdry) / eq.rdim;
		z_max = max(eq.zbdry) / eq.rdim;

		ChebyshevExpansion ch(n, psi, r_min, r_max, z_min, z_max);

		double d_psi = (eq.sibdry - eq.simagx) / (eq.nx - 1);

		for(size_t i = 0; i<N; i++){
			double r = r_min + (r_max - r_min) * i / (N - 1); // dimensionless
			for(size_t j = 0; j<N; j++){
				double z = z_min + (z_max - z_min) * j / (N - 1); // dimensionless

				// From \Psi definition
				Br(i, j) =  (sign ? -1.0 : 1.0) * ch.dy(r, z) / r;
				Bz(i, j) =  (sign ? 1.0 : -1.0) * ch.dx(r, z) / r;

				double psi_here = ch(r, z);
				size_t index = std::floor((psi_here - eq.simagx) / d_psi);
				// Keep index inside boundaries
				if (index <= 0) index = 1;
				if (index >= eq.nx - 1) index = eq.nx - 2;

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

class MagneticField{
	MagneticFieldMatrix& M;
	double B0;
public:
	MagneticField(MagneticFieldMatrix& B, double B_0) : M(B), B0(B_0) {}

	Vector3 operator()(Vector3 r, double /* t */ ){
		ScalarField MBr(M.Br, M.r_min, M.r_max, M.z_min, M.z_max);
		ScalarField MBt(M.Bt, M.r_min, M.r_max, M.z_min, M.z_max);
		ScalarField MBz(M.Bz, M.r_min, M.r_max, M.z_min, M.z_max);

		double x = r[0], y = r[2];
		double Br = six_point_formula(x, y, MBr) / B0;
		double Bt = six_point_formula(x, y, MBt) / B0;
		double Bz = six_point_formula(x, y, MBz) / B0;

		return Vector3 {Br, Bt, Bz};
	}
};

#endif // FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

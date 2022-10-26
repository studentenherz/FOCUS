#if !defined(FOCUS_INCLUDE_MAGNETIC_FIELD_HPP)
#define FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

#include "chebyshev.hpp"
#include "interpolations.hpp"
#include "types/equilibrium.hpp"
#include "types/matrix_2d.hpp"
#include "types/scalar_field.hpp"
#include "types/vector.hpp"

struct MagneticFieldMatrix{
	Matrix2D<double> psi;	// psi
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
	#ifdef __CUDACC__
	__host__
	#endif
	MagneticFieldMatrix(Equilibrium &eq, size_t n, size_t N, bool sign = true): psi(N, N), Br(N, N), Bt(N, N), Bz(N, N){

		// Dimensionless limits of the matrixes 
		double mr_min = eq.rleft / eq.rdim;
		double mr_max = eq.rleft / eq.rdim + 1;
		double mz_min = (eq.zmid - 0.5 * eq.zdim) / eq.rdim;
		double mz_max = (eq.zmid + 0.5 * eq.zdim) / eq.rdim;

		ScalarField raw_psi(eq.psi, mr_min, mr_max, mz_min, mz_max);
		
		// Dimensionless limits of plasma boundaries
		r_min = min(eq.rlim) / eq.rdim;
		r_max = max(eq.rlim) / eq.rdim;
		z_min = min(eq.zlim) / eq.rdim;
		z_max = max(eq.zlim) / eq.rdim;

		ChebyshevExpansion ch(n, raw_psi, r_min, r_max, z_min, z_max);

		for(size_t i = 0; i<N; i++){
			double r = r_min + (r_max - r_min) * i / (N - 1); // dimensionless r
			for(size_t j = 0; j<N; j++){
				double z = z_min + (z_max - z_min) * j / (N - 1); // dimensionless z

				// From \Psi definition
				Br(i, j) =  (sign ? -1.0 : 1.0) * (ch.dy(r, z) / (eq.bcentr * sqr(eq.rdim))) / r;
				Bz(i, j) =  (sign ? 1.0 : -1.0) * (ch.dx(r, z) / (eq.bcentr * sqr(eq.rdim))) / r;

				psi(i, j) = ch(r, z);
				double F = lagrange_interpolation_3(psi(i, j), eq.fpol, eq.simagx, eq.sibdry);

				// From F definition
				Bt(i, j) = (F / (eq.bcentr * eq.rdim)) / r;
			}
		}
	}

	/**
	 * Copy contructor for passing to the device
	 */
	#ifdef __CUDACC__
	__host__
	MagneticFieldMatrix(MagneticFieldMatrix& other): r_min(other.r_min), r_max(other.r_max), z_min(other.z_min), z_max(other.z_max)  {
		psi.construct_in_host_for_device(other.psi);
		Br.construct_in_host_for_device(other.Br);
		Bt.construct_in_host_for_device(other.Bt);
		Bz.construct_in_host_for_device(other.Bz);
	}
	#endif
};

class MagneticFieldFromMatrix{
	MagneticFieldMatrix& M;
	double _B0;
public:

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	MagneticFieldFromMatrix(MagneticFieldMatrix& B, double B_0) : M(B), _B0(B_0) {}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Vector3 operator()(Vector3 r, double /* t */ ){
		ScalarField MBr(M.Br, M.r_min, M.r_max, M.z_min, M.z_max);
		ScalarField MBt(M.Bt, M.r_min, M.r_max, M.z_min, M.z_max);
		ScalarField MBz(M.Bz, M.r_min, M.r_max, M.z_min, M.z_max);

		double x = r[0], y = r[2];
		double Br = six_point_formula(x, y, MBr);
		double Bt = six_point_formula(x, y, MBt);
		double Bz = six_point_formula(x, y, MBz);

		// This can't be done in the __device__
		// if (std::isnan(Br))
		// 	std::cerr << "Nan value of Br for r = " << r << '\n';
		// if (std::isnan(Bt))
		// 	std::cerr << "Nan value of Bt for r = " << r << '\n';
		// if (std::isnan(Bz))
		// 	std::cerr << "Nan value of Bz for r = " << r << '\n';

		return Vector3 {Br, Bt, Bz};
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double psi(Vector3 r, double /* t */ ){
		ScalarField Psi(M.psi, M.r_min, M.r_max, M.z_min, M.z_max);
		double x = r[0], y = r[2];
		return six_point_formula(x, y, Psi);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double B0(){
		return _B0;
	}
};

#ifndef __CUDACC__

class FineEquilibrium{
	Equilibrium& eq;
	ChebyshevExpansion ch;
	bool sign;
public:
	FineEquilibrium(Equilibrium& eq, size_t n, bool sign = true): 
	eq(eq), // Equilibrium 
	ch(			// Chebyshev Expansion
		n,		// order 
		ScalarField( // Scalar Field to be expanded
			eq.psi,			// matrix 
	 		eq.rleft / eq.rdim, // rmin
	 		eq.rleft / eq.rdim + 1, // rmax
	 		(eq.zmid - 0.5 * eq.zdim) / eq.rdim, // zmin
	 		(eq.zmid + 0.5 * eq.zdim) / eq.rdim), // xmax
		min(eq.rlim) / eq.rdim, // expansion rmin
		max(eq.rlim) / eq.rdim, // expansion rmax
		min(eq.zlim) / eq.rdim, // expansion zmin
		max(eq.zlim) / eq.rdim), // expansion zmax
		sign(sign)
	{}

	double Psi(double r, double z){
		return ch(r, z);
	}

	double F(double r, double z){
		double psi_here = Psi(r, z);
		return lagrange_interpolation_3(psi_here, eq.fpol, eq.simagx, eq.sibdry);
	}
	double Br(double r, double z){
		return (sign ? -1.0 : 1.0) * (ch.dy(r, z) / (eq.bcentr * sqr(eq.rdim))) / r;
	}

	double Bz(double r, double z){
		return (sign ? 1.0 : -1.0) * (ch.dx(r, z) / (eq.bcentr * sqr(eq.rdim))) / r;
	}


	double Bt(double r, double z){
		return (F(r, z) / (eq.bcentr * sqr(eq.rdim))) / r;
	}

	Vector3 B(double r, double z){
		return {Br(r, z), Bt(r, z), Bz(r, z)};
	}

	double B0(){
		return eq.bcentr;
	}
};

class MagneticField{
	FineEquilibrium& fineq;
public:
	MagneticField(FineEquilibrium& eq): fineq(eq) {}

	Vector3 operator()(Vector3 r, double /* t */){
		return fineq.B(r[0], r[2]);
	}

	double B0(){
		return fineq.B0();
	}
};

#endif // __CUDACC__

#endif // FOCUS_INCLUDE_MAGNETIC_FIELD_HPP

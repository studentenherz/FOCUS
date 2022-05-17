#if !defined(FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP)
#define FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP

#include "types/array.hpp"
#include "types/matrix_2d.hpp"

/** 
 * Struct to store the equilibrium data from
 * EFIT in the G-EQDSK files
 */
struct Equilibrium{
	int idnum;						// Id of file?
	size_t nx, ny; 				// Number of points in the R-z grid
	double rdim, zdim;		// Dimesions represented by the grid [meters]
	double rcentr;				// Reference value for R
	double bcentr;				// Vacuum toroidal magnetic field at rcentr
	double rleft;					// R at left (inner) boundary
  double zmid;					// z at middle of domain
  double rmagx, zmagx;	// R,z at magnetic axis (O-point)
  double simagx;        // Poloidal flux psi at magnetic axis
  double sibdry;        // Poloidal flux psi at plasma boundary
  double cpasma;        // Plasma current [Amperes]


	Array<double> fpol;		// 1D array of f(psi)=R*Bt  [meter-Tesla]
	Array<double> pres;		// 1D array of p(psi) [Pascals]
	Array<double> qpsi;		// 1D array of q(psi)

	Matrix2D<double> psi;	// 2D array (nx,ny) of poloidal flux

	// Boundaries description, optional
	size_t nbdry, nlim;					// Number of points of plasma and wall boundaries
	Array<double> rbdry, zbdry;	// 1D array of q(psi)
	Array<double> rlim, zlim;		// 1D array of q(psi)


	/**
	 * Default constructor 
	 * @param id idnum
	 * @param nx Number of points in R
	 * @param ny Number of points in z
	 */
	Equilibrium(int id = 0, size_t nx = 0, size_t ny = 0): idnum(id), nx(nx), ny(ny), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny), nbdry(0), nlim(0) {}
};

#endif // FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP



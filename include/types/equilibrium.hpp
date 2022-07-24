#if !defined(FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP)
#define FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP

#include "types/array.hpp"
#include "types/matrix_2d.hpp"

/** 
 * Struct to store the equilibrium data from
 * EFIT in the G-EQDSK files
 */
struct Equilibrium{
	/**
	 * @name Single parameters
	 */
	///@{
	int idnum;						///< Id of file
	size_t nx;						///< Number of points in the R direction
	size_t ny; 						///< Number of points in the z direction
	double rdim;					///< Dimensions of the represented grid in R (meters)
	double zdim;					///< Dimensions of the represented grid in z (meters)
	double rcentr;				///< Reference value for R (meters)
	double bcentr;				///< Vacuum toroidal magnetic field at `rcentr` (Tesla)
	double rleft;					///< R at left (inner) boundary (meters)
  double zmid;					///< z at middle of domain (meters)
  double rmagx;					///< R at magnetic axis (O-point) (meters)
	double zmagx;					///< z at magnetic axis (O-point) (meters)
  double simagx;        ///< Poloidal flux `psi` at magnetic axis
  double sibdry;        ///< Poloidal flux `psi` at plasma boundary
  double cpasma;        ///< Plasma current (Amperes)
	///@}

	/**
	 * @name Profiles
	 * One dimensional profiles as function of `psi` where the first value corresponds to `simagx` and the last to `sibdry`.
	 */
	///@{
	Array<double> fpol;		///< `F(psi) = R Bt` (meter-Tesla).
	Array<double> pres;		///< p(psi) (Pascal)
	Array<double> qpsi;		///< q(psi)
	///@}

	/**
	 * @name Poloidal flux
	 */
	///@{
	Matrix2D<double> psi;	///< 2D array (nx,ny) of poloidal flux
	///@}

	/**
	 * @name Boundaries description (optional)
	 */
	///@{
	size_t nbdry;								///< Number of points of plasma boundaries
	size_t nlim;								///< Number of points of wall limits
	Array<double> rbdry;				///< Plasma boundary R (meters)
	Array<double> zbdry;				///< Plasma boundary z (meters)
	Array<double> rlim;				///< Wall limits R (meters)
	Array<double> zlim;				///< Wall limits z (meters)
	///@}


	/**
	 * Default constructor 
	 * @param id idnum
	 * @param nx Number of points in R
	 * @param ny Number of points in z
	 */
	Equilibrium(int id = 0, size_t nx = 0, size_t ny = 0): idnum(id), nx(nx), ny(ny), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny), nbdry(0), nlim(0) {}
};

#endif // FOCUS_INCLUDE_TYPES_EQUILIBRIUM_HPP



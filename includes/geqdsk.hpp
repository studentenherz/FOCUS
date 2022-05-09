/**
 * Read geqdsk (G formatted EQuilibrium DiSK) files
 * 
 * This replicates what's implemented here:
 * https://github.com/bendudson/freegs/blob/master/freegs/_geqdsk.py
 * and here:
 * https://github.com/bendudson/pyTokamak/blob/master/tokamak/formats/geqdsk.py
 */
#if !defined(FOCUS_INCLUDES_GEQDSK_HPP)
#define FOCUS_INCLUDES_GEQDSK_HPP

#include <regex>
#include <fstream>
#include <string>
#include <iostream>

#include "types/array.hpp"
#include "types/matrix2D.hpp"

/** 
 * Struct to store the equilibrium data from
 * the G-EQDSK file
 */
struct Equilibrium{
	int idnum;						///< Id of file?
	size_t nx, ny; 				///< Number of points in the R-z grid
	double rdim, zdim;		///< Dimesions represented by the grid [meters]
	double rcentr;				///< Reference value for R
	double bcentr;				///< Vacuum toroidal magnetic field at rcentr
	double rleft;					///< R at left (inner) boundary
  double zmid;					///< z at middle of domain
  double rmagx, zmagx;	///< R,z at magnetic axis (O-point)
  double simagx;        ///< Poloidal flux psi at magnetic axis
  double sibdry;        ///< Poloidal flux psi at plasma boundary
  double cpasma;        ///< Plasma current [Amperes]

	Array<double> fpol;		///< 1D array of f(psi)=R*Bt  [meter-Tesla]
	Array<double> pres;		///< 1D array of p(psi) [Pascals]
	Array<double> qpsi;		///< 1D array of q(psi)

	Matrix2D<double> psi;	///< 2D array (nx,ny) of poloidal flux

	Equilibrium(int id, size_t nx, size_t ny): idnum(id), nx(nx), ny(nx), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny) {}
};

void read_eqdsk(const char *filename){

	std::string s = "TRXPL 27Aug2021 D3D D3D.13 153072G64 t~   3.4000   0 109 109  0.148607998E+01 0.294408007E+01";
	std::regex rexp("[+-]?\\d*[\\.]?\\d+(?:[Ee][+-]?\\d+)?");
	std::smatch m;

	 while (std::regex_search (s,m,rexp)) {
    std::cout << m[0] << "\n";
    s = m.suffix().str();
  }

}

#endif // FOCUS_INCLUDES_GEQDSK_HPP

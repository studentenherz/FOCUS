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
#include <sstream>

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

	/**
	 * Default constructor 
	 * @param id idnum
	 * @param nx Number of points in R
	 * @param ny Number of points in z
	 */
	Equilibrium(int id = 0, size_t nx = 0, size_t ny = 0): idnum(id), nx(nx), ny(ny), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny) {}
};

/**
 * Class to tokenize a stream based on a regex pattern
 */
template <typename stream_type>
class Tokenizer{
	std::regex e;
	std::string line;
	std::smatch m;
public:
	/**
	 * @param rexp Regular expresion
	 */
	Tokenizer(const char rexp[]): e(rexp) {}

	/**
	 * Get next token from stream
	 * @param is stream to tokenize
	 * @param[out] token string to put the token in
	 * @return true if new token, false otherwise
	 */
	bool next(stream_type& is, std::string& token){
		if (!std::regex_search (line,m,e)){
			if (!std::getline(is, line))
				return false;
			if (!std::regex_search (line,m,e))
				return false;
		}
    token = m[0];
		line = m.suffix().str();
		return true;
	}
};

/**
 * Read equilibrium G-EQDSK file
 * @param filename name of file
 * @return Equilibrium with the read data
 */
Equilibrium read_eqdsk(const char *filename){
	std::ifstream fi(filename);
	if(!fi.is_open()){
		std::cerr << "Couldn't open file " << filename << '\n';
		return Equilibrium(-1, 0, 0);
	}

	// First line contains a description and the last three tokens are
	// ... idnum nx ny
	std::string fline;
	std::getline(fi, fline);
	
	size_t i = 0;
	std::string tokens[3];
	std::istringstream iss(fline);
	while(iss >> tokens[i % 3]) i++;

	if (i < 3){
		std::cerr << "Error in" << filename << ". Expected at least 3 values in the first line.\n";
		return Equilibrium(-1, 0, 0);
	}

	int idnum = std::stoi(tokens[i % 3]);
	size_t nx = std::stoul(tokens[(i + 1) % 3]);
	size_t ny = std::stoul(tokens[(i + 2) % 3]);

	Equilibrium eq(idnum, nx, ny);

	eq.bcentr = 1.123;
	eq.fpol[12] = 23;
	eq.psi(45, 67) = 312.098;

	return eq;
}

#endif // FOCUS_INCLUDES_GEQDSK_HPP

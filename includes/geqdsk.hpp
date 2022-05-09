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

	Equilibrium(int id, size_t nx, size_t ny): idnum(id), nx(nx), ny(nx), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny) {}
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

void read_eqdsk(const char *filename){


	std::ifstream fi(filename);
	std::string token;
	std::getline(fi, token);
	std::cout << token << '\n';


	std::cout << "\nNow from Tokenizer\n";
	Tokenizer<std::ifstream> tk("[+-]?\\d*[\\.]?\\d+(?:[Ee][+-]?\\d+)?"); // captures any number;
	while(tk.next(fi, token)) std::cout << token << '\n';
}

#endif // FOCUS_INCLUDES_GEQDSK_HPP

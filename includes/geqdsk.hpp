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

	// Boundaries description, optional
	size_t nbdry, nlim;					///< Number of points of plasma and wall boundaries
	Array<double> rbdry, zbdry;	///< 1D array of q(psi)
	Array<double> rlim, zlim;		///< 1D array of q(psi)


	/**
	 * Default constructor 
	 * @param id idnum
	 * @param nx Number of points in R
	 * @param ny Number of points in z
	 */
	Equilibrium(int id = 0, size_t nx = 0, size_t ny = 0): idnum(id), nx(nx), ny(ny), fpol(nx), pres(nx), qpsi(nx), psi(nx, ny), nbdry(0), nlim(0) {}
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

	// First line of the file contains a description and the last three tokens are
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

	// Equilibrium variable to hold the data
	Equilibrium eq(idnum, nx, ny);

	// Give the 
	Tokenizer<std::ifstream> tk("[+-]?\\d*[\\.]?\\d+(?:[Ee][+-]?\\d+)?"); // captures any number;
	std::string token;

	// Next four lines of the file contain the experiment and tokamak characteristics
	if (tk.next(fi, token)) eq.rdim = std::stod(token);
	if (tk.next(fi, token)) eq.zdim = std::stod(token);
	if (tk.next(fi, token)) eq.rcentr = std::stod(token);
	if (tk.next(fi, token)) eq.rleft = std::stod(token);
	if (tk.next(fi, token)) eq.zmid = std::stod(token);

	if (tk.next(fi, token)) eq.rmagx = std::stod(token);
	if (tk.next(fi, token)) eq.zmagx = std::stod(token);
	if (tk.next(fi, token)) eq.simagx = std::stod(token);
	if (tk.next(fi, token)) eq.sibdry = std::stod(token);
	if (tk.next(fi, token)) eq.bcentr = std::stod(token);

	if (tk.next(fi, token)) eq.cpasma = std::stod(token);
	if (tk.next(fi, token)) eq.simagx = std::stod(token);
	if (tk.next(fi, token)); // here lies a dumb value
	if (tk.next(fi, token)) eq.rmagx = std::stod(token);
	if (tk.next(fi, token)) // here lies a dumb value

	if (tk.next(fi, token)) eq.zmagx = std::stod(token);
	if (tk.next(fi, token)); // here lies a dumb value
	if (tk.next(fi, token)) eq.sibdry = std::stod(token);
	if (tk.next(fi, token)); // here lies a dumb value
	if (tk.next(fi, token)); // here lies a dumb value

	// Read arrays
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) eq.fpol[i] = std::stod(token);
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) eq.pres[i] = std::stod(token);
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)); // no idea what are this values
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)); // no idea what are this values

	// Read matrix in "the natural" order
	for (size_t j = 0; j < eq.ny; j++)
		for (size_t i = 0; i < eq.nx; i++)
			if(tk.next(fi, token)) eq.psi(i, j) = std::stod(token);


	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) eq.qpsi[i] = std::stod(token);

	// Boundary and limits
	if (tk.next(fi, token)) eq.nbdry = std::stoul(token);
	if (tk.next(fi, token)) eq.nlim = std::stoul(token);

	if (eq.nbdry > 0) {
		eq.rbdry.resize(eq.nbdry);
		eq.zbdry.resize(eq.nbdry);
	
		for (size_t i = 0; i < eq.nbdry; i++){
			if (tk.next(fi, token)) eq.rbdry[i] = std::stod(token);
			if (tk.next(fi, token)) eq.zbdry[i] = std::stod(token);
		}
	}

	if (eq.nlim > 0) {
		eq.rlim.resize(eq.nlim);
		eq.zlim.resize(eq.nlim);
	
		for (size_t i = 0; i < eq.nlim; i++){
			if (tk.next(fi, token)) eq.rlim[i] = std::stod(token);
			if (tk.next(fi, token)) eq.zlim[i] = std::stod(token);
		}
	}

	return eq;
}

#endif // FOCUS_INCLUDES_GEQDSK_HPP

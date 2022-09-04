/**
 * Read geqdsk (G formatted EQuilibrium DiSK) files
 * 
 * This replicates what's implemented here:
 * https://github.com/bendudson/freegs/blob/master/freegs/_geqdsk.py
 * and here:
 * https://github.com/bendudson/pyTokamak/blob/master/tokamak/formats/geqdsk.py
 */
#if !defined(FOCUS_INCLUDE_FORMATS_GEQDSK_HPP)
#define FOCUS_INCLUDE_FORMATS_GEQDSK_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

#include "types/equilibrium.hpp"
#include "formats/regex_tokenizer.hpp"


/**
 * Read equilibrium G-EQDSK file
 * @param filename name of file
 * @return Equilibrium with the read data
 */
Equilibrium read_geqdsk(std::string filename){
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
	if (tk.next(fi, token)) {} // here lies a dumb value
	if (tk.next(fi, token)) eq.rmagx = std::stod(token);
	if (tk.next(fi, token)) {} // here lies a dumb value

	if (tk.next(fi, token)) eq.zmagx = std::stod(token);
	if (tk.next(fi, token)) {} // here lies a dumb value
	if (tk.next(fi, token)) eq.sibdry = std::stod(token);
	if (tk.next(fi, token)) {} // here lies a dumb value
	if (tk.next(fi, token)) {} // here lies a dumb value

	// Read arrays
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) eq.fpol[i] = std::stod(token);
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) eq.pres[i] = std::stod(token);
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) {} // no idea what are this values
	for (size_t i = 0; i < eq.nx; i++) if(tk.next(fi, token)) {} // no idea what are this values

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

#endif // FOCUS_INCLUDE_FORMATS_GEQDSK_HPP

#include <fstream>
#include <string>
#include <cassert>

#include "formats/geqdsk.hpp"
#include "formats/matrix.hpp"
#include "cxxopts.hpp"
#include "magnetic_field.hpp"
#include "types/vector.hpp"

__global__ 
void kernel(Equilibrium eq, MagneticFieldMatrix B_matrix, double *psi_sum, double *fpol_sum, int *idnum, Vector3* sum_b){
	*idnum = eq.idnum;

	*fpol_sum = 0;
	for (size_t i = 0; i < eq.fpol.size(); i++)
		*fpol_sum += eq.fpol[i];

	*psi_sum = 0;
	for(size_t i = 0; i < eq.psi.shape().first; i++)
		for(size_t j = 0; j < eq.psi.shape().second; j++)
			*psi_sum += eq.psi(i, j);

	MagneticFieldFromMatrix B(B_matrix, eq.bcentr);
	*sum_b = {0, 0, 0};
		for(size_t i = 0; i < 100; i++)
			for(size_t j = 0; j <  100; j++){
				double r = B_matrix.r_min + i * (B_matrix.r_max - B_matrix.r_min) / 100;
				double z = B_matrix.z_min + j * (B_matrix.z_max - B_matrix.z_min) / 100;
				*sum_b = *sum_b + B({r, 0, z}, 0);
			}
}

int main(int argc, char* argv[]){
	
	cxxopts::options options("geqdsk", "Test geqdsk read");

	options.add_options()
		("file", "", cxxopts::value<std::string>())
		("h,help", "Show this help message");

	options.positional_help("<G-EQDSK input file>");
	options.parse_positional({"file"});

	try{
		auto result = options.parse(argc, argv);

		if (result.count("help")){
			std::cout << options.help() << std::endl;
			return 0;
		}

		Equilibrium eq = read_geqdsk(result["file"].as<std::string>());
		std::cout << eq.idnum << '\n';
		std::cout << eq.nx << '\n';
		std::cout << eq.ny << '\n';
		std::cout << eq.rdim << '\n';
		std::cout << eq.zdim << '\n';
		std::cout << eq.rcentr << '\n';
		std::cout << eq.rleft << '\n';
		std::cout << eq.zmid << '\n';
		std::cout << eq.rmagx << '\n';
		std::cout << eq.zmagx << '\n';
		std::cout << eq.simagx << '\n';
		std::cout << eq.sibdry << '\n';
		std::cout << eq.bcentr << '\n';
		std::cout << eq.cpasma << '\n';
		std::cout << eq.simagx << '\n';
		std::cout << eq.rmagx << '\n';
		std::cout << eq.zmagx << '\n';
		std::cout << eq.sibdry << '\n';

		dump("psi_from_geqdsk.dat", eq.psi, false);
		
		std::ofstream fo1("psi.dat");
		for (size_t j = 0; j < eq.ny; j++)
			for (size_t i = 0; i < eq.nx; i++){
				double r = eq.rleft + i * eq.rdim / eq.nx;
				double z = eq.zmid - eq.zdim / 2 + j * eq.zdim / eq.ny;
				fo1 << r << ' ' << z << ' ' << eq.psi(i, j) << '\n';
			}


		// Output limit
		std::ofstream fo("limit.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			fo << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';

		std::ofstream fo2("boundary.dat");
		for (size_t i = 0; i < eq.nlim; i++)
			fo2 << eq.rlim[i] << ' ' << eq.zlim[i] << '\n';
		
		double fpol_sum = 0;
		for (size_t i = 0; i < eq.fpol.size(); i++)
			fpol_sum += eq.fpol[i];

		double psi_sum = 0;	
		for(size_t i = 0; i < eq.psi.shape().first; i++)
			for(size_t j = 0; j < eq.psi.shape().second; j++)
				psi_sum += eq.psi(i, j);

		// Equilibrium dEq;
		// dEq.construct_in_host_for_device(eq);

		MagneticFieldMatrix B_matrix(eq, 26, 400);
		MagneticFieldFromMatrix B(B_matrix, eq.bcentr);

		Vector3 sum_b = {0, 0, 0};
		for(size_t i = 0; i < 100; i++)
			for(size_t j = 0; j <  100; j++){
				double r = B_matrix.r_min + i * (B_matrix.r_max - B_matrix.r_min) / 100;
				double z = B_matrix.z_min + j * (B_matrix.z_max - B_matrix.z_min) / 100;
				sum_b = sum_b + B({r, 0, z}, 0);
			}

		double *d_psi_sum, *d_fpol_sum;
		int *d_idnum;
		Vector3 *d_sum_b;

		cudaMalloc(&d_psi_sum, sizeof(double));
		cudaMalloc(&d_fpol_sum, sizeof(double));
		cudaMalloc(&d_idnum, sizeof(int));
		cudaMalloc(&d_sum_b, sizeof(Vector3));

		kernel<<<1, 1>>>(eq, B_matrix, d_psi_sum, d_fpol_sum, d_idnum, d_sum_b);

		double h_psi_sum, h_fpol_sum;
		int h_idnum;
		Vector3 h_sum_b;

		cudaMemcpy(&h_psi_sum, d_psi_sum, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_fpol_sum, d_fpol_sum, sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_idnum, d_idnum, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&h_sum_b, d_sum_b, sizeof(Vector3), cudaMemcpyDeviceToHost);

		assert(fpol_sum == h_fpol_sum);
		assert(psi_sum == h_psi_sum);
		assert(h_idnum == eq.idnum);
		std::cout << h_sum_b << '\t' << sum_b << '\n';

		return 0;
	}catch(cxxopts::option_error const& e){
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
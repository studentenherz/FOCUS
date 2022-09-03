#include <iostream>
#include <string>
#include <fstream>

#include "cxxopts.hpp"
#include "random.hpp"
#include "types/array.hpp"

__global__
void kernel_init_rand(PhiloxCuRand philox, unsigned long long seed){
  philox.init(seed);
}

__global__
void kernel_generate(PhiloxCuRand philox, Array<double> unif, Array<double> norm){
  size_t idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx >= unif.size()) return;
  unif[idx] = philox.uniform();
  norm[idx] = philox.normal();
}

// This is that need to be actually well randomized: the random nombers for each particle
// that will be in a single thread
__global__
void kernel_generate_serial(PhiloxCuRand philox, Array<double> unif, Array<double> norm){
  for (size_t i = 0; i<unif.size(); i++){
    unif[i] = philox.uniform();
    norm[i] = philox.normal();
  }
}

int main(int argc, char* argv[]){
  cxxopts::options options("random", "Test CUDA Philox random");

  options.add_options()
    ("n,number", "Number of points to calculate", cxxopts::value<size_t>()->default_value("1000"))
    ("N,normal", "Normal destribution out file name", cxxopts::value<std::string>()->default_value("norm.dat"))
    ("U,uniform", "Uniform destribution out file name", cxxopts::value<std::string>()->default_value("unif.dat"))
    ("s,serial", "Calculate the distribution in same thread")
    ("h,help", "Show this help message");

  try{
    auto result = options.parse(argc, argv);
    
    if (result.count("help")){
      std::cout << options.help() << std::endl;
      return EXIT_SUCCESS;
    }

    size_t Ntotal = result["number"].as<size_t>();
    size_t threadsPerBlock = 1000;
    size_t nBlocks = Ntotal / threadsPerBlock;

    if(nBlocks * threadsPerBlock < Ntotal) nBlocks++;

    size_t Nthreads = nBlocks * threadsPerBlock;

    PhiloxCuRand philox(Nthreads);
    kernel_init_rand<<<nBlocks, threadsPerBlock>>>(philox, 1234);

    Array<double> h_unif(Ntotal);
    Array<double> h_norm(Ntotal);

    Array<double> d_unif, d_norm;
    d_unif.construct_in_host_for_device(h_unif);
    d_norm.construct_in_host_for_device(h_norm);

    if (result.count("serial"))
      kernel_generate_serial<<<1, 1>>>(philox, d_unif, d_norm);
    else
      kernel_generate<<<nBlocks, threadsPerBlock>>>(philox, d_unif, d_norm);

    h_unif.copy_to_host_from_device(d_unif);
    h_norm.copy_to_host_from_device(d_norm);

    std::ofstream funif(result["uniform"].as<std::string>());
    std::ofstream fnorm(result["normal"].as<std::string>());

    for(size_t i =0 ; i < Ntotal; i++){
      funif << h_unif[i] << '\n';
      fnorm << h_norm[i] << '\n';
    }

    funif.close();
    fnorm.close();

    return EXIT_SUCCESS;
  }catch(cxxopts::parse_error const& e){
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

}
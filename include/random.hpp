#if !defined(FOCUS_INCLUDE_RANDOM_HPP)
#define FOCUS_INCLUDE_RANDOM_HPP

#include <cmath>
#include <curand_kernel.h>

#include "util.hpp"

/**
 * Ran2 implementation from Numerical Recipes E.3
 */
class Ran2
{
	unsigned long long u, v, w;
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Ran2(unsigned long long j) : v(4101842887655102017LL), w(1){
		u = j ^ v;
		int64();
		v = u;
		int64();
		w = v;
		int64();
	}

	/**
	 * Get random unsigned int 64 bits (uniform distribution)
	 * @return random unit64
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	inline unsigned long long int64(){
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17;
		v ^= v << 31;
		v ^= v >> 8;
		w = 4294957665U * (w & 0xffffffff) + (w >> 32);
		unsigned long long x = u ^ (u << 21);
		x ^= x >> 35;
		x ^= x << 4;
		return (x + v) ^ w;
	}

	/**
	 * Get random unsigned int 32 bits (uniform distribution)
	 * @return random unit32
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	inline unsigned int int32() { return (unsigned int)int64(); }
	
	/**
	 * Get random float 64 bits (uniform distribution)
	 * @return random float64
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	
	/**
	 * Get random float 64 bits in a range (uniform distribution)
	 * @param xmin lower limit
	 * @param xmax higher limit
	 * @return random float64 between xmin and xmax
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double random(double xmin, double xmax){
		return xmin + doub() * (xmax - xmin);
	}
};

/**
 * Normal distribution random number generator
 */
class NormalRand{
	double sigma;	// std deviation
	double mu;		// mean
	unsigned long long seed;
	Ran2 ran;			// uniform random generator
public:
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	NormalRand(unsigned long long seed_, double sigma_ = 1.0, double mu_ = 0): sigma(sigma_), mu(mu_), seed(seed_), ran(seed) {};
	
	/**
	 * Get random float 64 from normal distribution
	 * @return random float64
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double operator()(){
		// This uses Box-Muller algorithm
		double x1, x2;
		x1 = ran.doub(); // 0 to 1
		x2 = ran.doub(); // 0 to 1
		return sigma * sqrt(-2.0 * log(x1)) * cos(two_pi * x2) + mu;
	}
};

#ifdef __CUDACC__

class PhiloxCuRand{
	curandStatePhilox4_32_10_t *d_states;
	size_t _n;
public:
	__host__
	PhiloxCuRand(size_t n): _n(n) {
		cudaMalloc(&d_states, _n * sizeof(curandStatePhilox4_32_10_t));
	}

	__device__
	void init(unsigned long long seed, unsigned long long subsequence = 0LL, unsigned long long offset = 0LL){
		size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
		curand_init(seed, idx, offset, &d_states[idx]);
	}

	__device__
	double uniform(){
		size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
		return curand_uniform_double(&d_states[idx]);
	}

	__device__
	double normal(){
		size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
		return curand_normal_double(&d_states[idx]);
	}
};

#endif // __CUDACC__

#endif // FOCUS_INCLUDE_RANDOM_HPP

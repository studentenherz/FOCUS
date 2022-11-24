#if !defined(FOCUS_INCLUDE_TYPES_PLASMA_HPP)
#define FOCUS_INCLUDE_TYPES_PLASMA_HPP

#include <vector>
#include <string>

#include "interpolations.hpp"
#include "types/array.hpp"
#include "types/matrix_2d.hpp"
#include "particle.hpp"
#include "util.hpp"

struct Plasma{
	int shot;
	size_t nexp;
	size_t nion;
	double masse;
	double ze;
	double logl_prefactor;

	Array<double> mass;
	Array<double> z;

	Array<double> polflux;

	Array<double> ne;
	Matrix2D<double> ni;	

	Array<double> te;
	Matrix2D<double> ti;

	Plasma(int shot, size_t nexp, size_t nion): shot(shot), nexp(nexp), nion(nion), mass(nion), z(nion), polflux(nexp), ne(nexp), ni(nion, nexp), te(nexp), ti(nion, nexp) {}

	#ifdef __CUDACC__
	__host__
	void construct_in_host_for_device(Plasma& other){
		logl_prefactor = other.logl_prefactor;
		shot = other.shot;
		nexp = other.nexp;
		nion = other.nion;
		masse = other.masse;
		ze = other.ze;

		mass.construct_in_host_for_device(other.mass);
		z.construct_in_host_for_device(other.z);

		polflux.construct_in_host_for_device(other.polflux);

		ne.construct_in_host_for_device(other.ne);
		te.construct_in_host_for_device(other.te);

		ni.construct_in_host_for_device(other.ni);
		ti.construct_in_host_for_device(other.ti);
	}
	#endif

	/**
	 * Interpolated temperature
	 * @param psi poloidal flux
	 * @param t time
	 * @return interpolated temperature
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double T(double psi, double /* t */){
		if (index < nion)
			return lagrange_interpolation_3(psi, polflux, ti, index);
		if (index == nion)
			return lagrange_interpolation_3(psi, polflux, te);
		return 0;
	}

	/**
	 * Interpolated density
	 * @param psi poloidal flux
	 * @param t time
	 * @return interpolated density
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double n(double psi, double /* t */){
		if (index < nion)
			return lagrange_interpolation_3(psi, polflux, ni, index);
		if (index == nion)
			return lagrange_interpolation_3(psi, polflux, ne);
		return 0;
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double m(){
		if (index > nion) return 0;
		return (index < nion ? mass[index] : masse);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double q(){
		if (index > nion) return 0;
		return (index < nion ? z[index] : ze);
	}

	#ifdef __CUDACC__
	__host__ __device__
	#endif
	double logl(Particle a, double psi, double t){
		if (index > nion) return 0;
		double m = (index < nion ? mass[index] : masse);
		double q = (index < nion ? z[index] : ze);

		return logl_prefactor + log(T(psi, t) / (n(psi, t) * sqr(1 + m/a.m) * sqr(a.q) * sqr(sqr(q))));
	}

	/**
	 * Size of plasma
	 * @return number of particle species
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	size_t size(){
		return nion + 1;
	}

	/**
	 * This allows to do this:
	 * ```
	 * for(size_t i = 0; i < plasma.size(); i++)
	 * 		double t = plasma[i].T(psi, t);
	 * ```
	 *
	 * @param i index of particle of plasma
	 * @return plasma with T and n evaluable for i-th particle species
	 */
	#ifdef __CUDACC__
	__host__ __device__
	#endif
	Plasma operator[](size_t i){
		index = i;
		return *this;
	}

private:
	size_t index;
};

#endif // FOCUS_INCLUDE_TYPES_PLASMA_HPP

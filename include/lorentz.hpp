#if !defined(FOCUS_INCLUDE_LORENTZ_HPP)
#define FOCUS_INCLUDE_LORENTZ_HPP

#include "types/vector.hpp"

/**
 * A vector field class that returns zero
 */
class NullVectorField {
public:
	NullVectorField(){}
	Vector3 operator()(Vector3 /* r */, double /* t */) {Vector3 zero; return zero;}
} null_vector_field; // Null Vector Field object

/**
 * A force class that returns zero
 */
class NullForce{
public:
	NullForce(){}
	Vector3 operator()(State /* x */, double /* t */) {Vector3 zero; return zero;}
} null_force; // Null Force object


/**
 * Motion equations with Lorentz force applied
 */
template<typename force_type, typename magnetic_field_type, typename electric_field_type = magnetic_field_type>
class MotionEquation{
	const double gam;			// dimensionless factor
	magnetic_field_type B;	// magnetic induction field
	electric_field_type E;	// electric field
	force_type F;					// other forces
public:
	/**
	 * Constructor of Motion Equation
	 * @param _gam dimensionless gamma factor
	 * @param _B magnetic field, callable B(Vector3 r, double t)
	 * @param _E electric field, callable E(Vector3 r, double t)
	 * @param _F additional force, callable F(State x, double t)
	 */
	MotionEquation(double _gam, magnetic_field_type _B = null_vector_field, electric_field_type _E = null_vector_field, force_type _F = null_force): gam(_gam), B(_B), E(_E), F(_F) {}
	
	/**
	 * Equation system for integration using the Lorentz force
	 * @param x current state
	 * @param dxdt state derivative
	 * @param t current time
	 */
	void operator()(const State &x, State &dxdt, double t ){

		Vector3 r = get_position(x);
		Vector3 b = B(r, t);
		Vector3 e = E(r, t);
		Vector3 f = F(x, t);	// other forces

		dxdt[0] = gam * x[3];															// d(rho)/dt = v_rho
		dxdt[1] = gam * x[4] / x[0];											// d(theta)/dt = v_theta / rho
		dxdt[2] = gam * x[5];															// dz/dt = v_z
		dxdt[3] = f[0] + x[4] * b[2] - x[5] * b[1] + e[0] + gam * x[4] * x[4] / x[0];		// v_rho
		dxdt[4] = f[1] + x[5] * b[0] - x[3] * b[2] + e[1] - gam * x[3] * x[4] / x[0];		// v_theta
		dxdt[5] = f[2] + x[3] * b[1] - x[4] * b[0] + e[2];															// v_z
	}
};

#endif // FOCUS_INCLUDE_LORENTZ_HPP

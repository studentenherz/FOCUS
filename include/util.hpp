#if !defined(FOCUS_INCLUDE_UTIL_HPP)
#define FOCUS_INCLUDE_UTIL_HPP

#include <cmath>

const double pi(3.14159265358979323846);
const double two_pi(6.283185307179586477);
const double two_over_sqrt_pi(1.1283791670955125738);

#ifdef __CUDACC__
__host__ __device__
#endif
inline double sqr(double x){return x * x;}

#endif // FOCUS_INCLUDE_UTIL_HPP

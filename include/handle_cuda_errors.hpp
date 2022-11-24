#if !defined(FOCUS_INCLUDE_HANDLE_CUDA_ERRORS_HPP)
#define FOCUS_INCLUDE_HANDLE_CUDA_ERRORS_HPP

#include <iostream>

#ifdef __CUDACC__

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
	if (code != cudaSuccess){
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

#define propagateCUDAErr(ans) { \
	cudaError_t err = (ans); \
	if (err != cudaSuccess) return err; \
}

inline void checkCUDAError(const char* msg, const char *file, int line, bool abort=true){
	cudaError_t code = cudaGetLastError();
	if (code != cudaSuccess){
		fprintf(stderr, "%s: %d %s %s %s %d\n", msg, code, cudaGetErrorName(code), cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#define checkKernelErr() checkCUDAError("Error in kernel launch", __FILE__, __LINE__)

#endif

#endif // FOCUS_INCLUDE_HANDLE_CUDA_ERRORS_HPP

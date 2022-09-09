#if !defined(FOCUS_INCLUDE_CUDA_PARALELL_REDUCTION_HPP)
#define FOCUS_INCLUDE_CUDA_PARALELL_REDUCTION_HPP

#include "types/array.hpp"

static const size_t blockSize = 1024;
static const size_t gridSize = 24;

template<typename T>
__global__ 
void kernel_sum_commutative_multi_block(Array<T> gArr, Array<T> gOut){
	size_t thIdx = threadIdx.x;
	size_t gthIdx = thIdx + blockIdx.x * blockSize;
	const size_t gridSize = blockSize * gridDim.x;

	T sum = 0;
	for (size_t i = gthIdx; i < gArr.size(); i+= gridSize)
		sum += gArr[i];

	__shared__ T shArr[blockSize];
	shArr[thIdx] = sum;

	__syncthreads();
	for (size_t size = blockSize / 2; size > 0; size /= 2){
		if(thIdx < size)
			shArr[thIdx] += shArr[thIdx + size];
		__syncthreads();
	}
	
	if (thIdx == 0)
		gOut[blockIdx.x] = shArr[0];
}

template<typename T>
T paralell_sum(Array<T>& d_Arr){
	Array<T> h_Out(gridSize);
	Array<T> d_Out;
	d_Out.construct_in_host_for_device(h_Out);

	kernel_sum_commutative_multi_block<<<gridSize, blockSize>>>(d_Arr, d_Out); // Partial result
	kernel_sum_commutative_multi_block<<<1, blockSize>>>(d_Out, d_Out); // Result in d_Out[0]

	T sum = d_Out.from_device_at(0);
	d_Out.deallocate_cuda();
	return sum;
}

#endif // FOCUS_INCLUDE_CUDA_PARALELL_REDUCTION_HPP

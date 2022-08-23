#include <iostream>

#include "types/array.hpp"



__global__ void k_sum_array(int dArr[], size_t* n, int* s){
	Array<int> arr(dArr, *n);
	
	*s = 0;
	for(size_t i = 0; i < arr.size(); i++)
		*s += arr[i];
}

int main(){
	size_t n = 10;
	size_t* dn;
	int hArr[n];
	for(size_t i = 0; i < n; i++)
		hArr[i] = i;
	int* dArr;

	int hsum;
	int* dsum;
	cudaMalloc(&dsum, sizeof(int));
	cudaMalloc(&dn, sizeof(size_t));
	cudaMalloc(&dArr, sizeof(int) * n);

	cudaMemcpy(dArr, hArr, sizeof(int) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(dn, &n, sizeof(size_t), cudaMemcpyHostToDevice);

	k_sum_array<<<1, 1>>>(dArr, dn, dsum);

	cudaMemcpy(&hsum, dsum, sizeof(int), cudaMemcpyDeviceToHost);

	cudaFree(dsum);
	cudaFree(dArr);
	cudaFree(dn);
	
	std::cout << hsum << '\n';

	return 0;
}
#include <iostream>

#include "types/array.hpp"

__device__ Array<int> get_arr(size_t n){
	Array<int> arr(n);
	for(size_t i = 0; i < arr.size(); i++)
		arr[i] = i + 1;
	return arr; // move 
}

__global__ void k_sum_array(Array<int> arr, int* s){
	// Array<int> arr(dArr, *n); // construct for using the pointer from host and then don't try to deallocate
	Array<int> arr2 = get_arr(arr.size()); // move constructor and dealocate when leaving scope
	Array<int> arr3 = arr2; // copy constructor and dealocate when leaving scope
	Array<int> arr4 = arr; // copy constructor from another constructed from device, dealocate the new allocated when leaving scope

	*s = 0;
	for(size_t i = 0; i < arr.size(); i++)
		*s += arr[i] + arr2[i] + arr3[i] + arr4[i];
}


int main(){
		size_t n = 10;
		int raw_arr[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
		int hsum;
		int* dsum;
		
		cudaMalloc(&dsum, sizeof(int));

		Array<int> arr(raw_arr, n);

		k_sum_array<<<1, 1>>>(arr, dsum);

		cudaMemcpy(&hsum, dsum, sizeof(int), cudaMemcpyDeviceToHost);
		std::cout << hsum << '\n';
	return 0;
}
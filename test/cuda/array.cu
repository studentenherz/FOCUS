#include <iostream>

#include "types/array.hpp"

__device__ Array<int> get_arr(size_t n){
	Array<int> arr(n);
	for(size_t i = 0; i < arr.size(); i++)
		arr[i] = i;
	return arr; // move 
}

__global__ void k_sum_array(int dArr[], size_t* n, int* s){
	Array<int> arr(dArr, *n); // construct for using the pointer from host and then don't try to deallocate
	Array<int> arr2 = get_arr(*n); // move constructor and dealocate when leaving scope
	Array<int> arr3 = arr2; // copy constructor and dealocate when leaving scope
	Array<int> arr4 = arr; // copy constructor from another constructed from device, dealocate the new allocated when leaving scope

	*s = 0;
	for(size_t i = 0; i < arr.size(); i++)
		*s += arr[i] + arr2[i] + arr3[i] + arr4[i];
}

class DoTheWork{
	size_t hn;
	size_t* dn;
	int* hArr;
	int* dArr;
	int hsum;
	int* dsum;
public:
	DoTheWork(size_t n): hn(n) {
		cudaMalloc(&dsum, sizeof(int));
		cudaMalloc(&dn, sizeof(size_t));
		cudaMalloc(&dArr, sizeof(int) * n);

		hArr = new int[n];
		for(size_t i = 0; i < n; i++)
			hArr[i] = i;
		cudaMemcpy(dArr, hArr, sizeof(int) * n, cudaMemcpyHostToDevice);
		cudaMemcpy(dn, &n, sizeof(size_t), cudaMemcpyHostToDevice);
	}

	~DoTheWork(){
		cudaFree(dsum);
		cudaFree(dArr);
		cudaFree(dn);
	}

	void now_do_it(){
		k_sum_array<<<1, 1>>>(dArr, dn, dsum);
		cudaMemcpy(&hsum, dsum, sizeof(int), cudaMemcpyDeviceToHost);
		std::cout << hsum << '\n';
	}
};

int main(){
	
	DoTheWork yo(11);
	yo.now_do_it();

	return 0;
}
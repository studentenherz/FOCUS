#include <iostream>

#include "types/array.hpp"

__global__ void k_sum_array(int dArr[], size_t* n, int* s){
	Array<int> arr(dArr, *n);
	
	*s = 0;
	for(size_t i = 0; i < arr.size(); i++)
		*s += arr[i];
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
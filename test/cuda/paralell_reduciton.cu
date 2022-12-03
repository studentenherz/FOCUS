#include <iostream>

#include "types/array.hpp"
#include "cuda/paralell_reduction.hpp"

int main(){
	size_t N = 1000000;
	Array<unsigned long long> h_Arr(N);
	unsigned long long sum = 0;
	for (size_t i = 0; i < N; i++) {
		h_Arr[i] = i;
		sum += i;
	}
	Array<unsigned long long> d_Arr;
	d_Arr.construct_in_host_for_device(h_Arr);

	std::cout << "paralell_sum: " << paralell_sum(d_Arr) << '\n';
	std::cout << "serial_sum: " << sum << '\n';

	return EXIT_SUCCESS;
}
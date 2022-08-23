#include <iostream>

#include "types/vector.hpp"

__global__ void kernel(Vector3* vec){
	Vector3 zero;
	Vector3 a = {1.2, 4.0, -9}, b;
	b = {0, 1, 3};

	*vec = a + b;
}

int main(){
	using namespace std;
	
	Vector3 zero;
	Vector3 a = {1.2, 4.0, -9}, b;
	b = {0, 1, 3};

	cout << "zero = " << zero << '\n';
	cout << "a = " << a << '\n';
	cout << "b = " << b << '\n';

	cout << "2 * a = " << 2 * a << '\n';
	cout << "a * (-2) = " << a * (-2) << '\n';
	cout << "a / 2 = " << a / 2 << '\n';

	cout << "a + b = " << a + b << '\n';
	cout << "a - b = " << a - b << '\n';
	cout << "dot(a, b) = " << dot(a, b) << '\n';
	cout << "cross(a, b) = " << cross(a, b) << '\n';
	cout << "mod(a) = " << mod(a) << '\n';

	Vector3 h_vec;
	Vector3* d_vec;

	cudaMalloc(&d_vec, sizeof(Vector3));

	kernel<<<1, 1>>>(d_vec);

	cudaMemcpy(&h_vec, d_vec, sizeof(Vector3), cudaMemcpyDeviceToHost);

	cout << h_vec << '\n';
	return 0;
}
#include <iostream>

#include "types/vector.hpp"

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

	return 0;
}
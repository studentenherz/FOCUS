#include <iostream>

#include "types/matrix_2d.hpp"

int main(){
	Matrix2D<int> m(2, 2);
	m(2, 2) = 1; // This should err-out because is out of bounds

	std::cout << "Hehe, I got here (evil bug's laugh)\n";

	return 0;
}
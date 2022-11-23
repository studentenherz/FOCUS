#include <iostream>

#include "types/array.hpp"

int main(){
	Array<int> a(2);
	a[2] = 1; // This should err-out because is out of bounds

	std::cout << "Hehe, I got here (evil bug's laugh)\n";

	return 0;
}
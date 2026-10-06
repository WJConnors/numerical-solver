#include <iostream>
#include <cmath>
#include <numbers>

double u(double x) {
	return std::sin(x * std::numbers::pi);
}

int main() {
    std::cout << "u(0) = " << u(0) << std::endl;
    std::cout << "u(0.5) = " << u(0.5) << std::endl;
    std::cout << "u(1) = " << u(1) << std::endl;
}

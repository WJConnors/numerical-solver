#include <iostream>
#include <numbers>
#include <cmath>

#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

//To be used for final verification, assuming u is unknown otherwise
//u(x,y) = sin(pi.x).sin(pi.y)
double exact_solution(double x, double y) {
	return std::sin(pi * x) * std::sin(pi * y);
}
//Laplacian
// f(x,y) = -2.pi^2.sin(pi.x).sin(pi.y)
double exact_laplacian(double x, double y) {
	return -2 * pi * pi * exact_solution(x, y);
}

int main() {
	std::cout << "Hello, poisson" << '\n';
}

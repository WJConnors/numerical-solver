#include <iostream>
#include <numbers>
#include <cmath>

#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

//Domain
//0 <= x <= 1
//< <= y <= 1
constexpr double min = 0.0;
constexpr double max = 1.0;
//n-1 for number of intervals
constexpr int ns[] = {17, 33, 65, 129};

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

	int n = ns[0];
	
	double d = (max - min) / (n -1);

	Grid2D grid = make_grid2d(n, n);

	for (int x = 0; x < n; x++) {
		for (int y = 0; y < n; y++) {
			grid(x, y) = exact_solution(
				min + x * d,
				min + y * d
			);
		}
	}	
}

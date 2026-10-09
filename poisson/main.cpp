#include <iostream>
#include <numbers>
#include <cmath>
#include <utility>

#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

constexpr double min = 0.0;
constexpr double max = 1.0;
constexpr int n = 17;

constexpr double tolerance = 1e-8;
constexpr int max_iterations = 1000000;

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
	double d = (max - min) / (n - 1);
	Grid2D f = make_grid2d(n, n);
	for (int y = 0; y < n; y++) {
		for (int x = 0; x < n; x++) {
			f(x, y) = exact_laplacian(
				min + x * d,
				min + y * d
			);
		}
	}

	Grid2D u_old = make_grid2d(n, n);
	Grid2D u_new = make_grid2d(n, n);
	//Jacobi iteration
	for (int i = 0; i < max_iterations; i++) {
		double max_change = 0.0;
		std::swap(u_old, u_new);
		for (int y = 1; y < n - 1; y++) {
			for (int x = 1; x < n - 1; x++) {
				u_new(x, y) = 
					(u_old(x + 1, y) +
					u_old(x - 1, y) +
					u_old(x, y + 1) +
					u_old(x, y - 1) -
					d * d * f(x, y))
					/ 4.0;

				double change = std::abs(u_new(x, y) - u_old(x, y));
				if (change > max_change) max_change = change;
			}
		}
		if (max_change < tolerance) {
			std::cout << "Converged after " << i + 1 << " iterations\n";
			break;
		}
	}

	double max_error = 0.0;
	for (int y = 1; y < n - 1; y++) {
		for (int x = 1; x < n - 1; x++) {
			double error = std::abs(
				u_new(x, y) -
				exact_solution(
					min + x * d,
					min + y * d
				));
			if (error > max_error) max_error = error;
		}
	}

	std::cout << max_error << '\n';
	
}

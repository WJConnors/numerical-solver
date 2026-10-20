#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include <utility>

#include "finite_difference.hpp"
#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

constexpr double min = 0.0;
constexpr double max = 1.0;
constexpr int ns[] = {17, 33, 65, 129};

constexpr int max_iterations = 1000000;
constexpr double residual_tolerance = 1e-8;

// To be used for final verification, assuming u is unknown otherwise
// u(x,y) = sin(pi.x).sin(pi.y)
double exact_solution(double x, double y) {
	return std::sin(pi * x) * std::sin(pi * y);
}
// Laplacian
//  f(x,y) = -2.pi^2.sin(pi.x).sin(pi.y)
double exact_laplacian(double x, double y) {
	return -2 * pi * pi * exact_solution(x, y);
}

int main() {
	for (int n : ns) {
		auto start = std::chrono::steady_clock::now();
		std::cout << "\nn = " << n << '\n';
		double d = (max - min) / (n - 1);
		Grid2D f = make_grid2d(n, n);
		for (int y = 0; y < n; y++) {
			for (int x = 0; x < n; x++) {
				f(x, y) = exact_laplacian(min + x * d, min + y * d);
			}
		}

		Grid2D u_old = make_grid2d(n, n);
		Grid2D u_new = make_grid2d(n, n);
		// Jacobi iteration
		for (int i = 0; i < max_iterations; i++) {
			std::swap(u_old, u_new);
			for (int y = 1; y < n - 1; y++) {
				for (int x = 1; x < n - 1; x++) {
					u_new(x, y) =
					    (u_old(x + 1, y) + u_old(x - 1, y) + u_old(x, y + 1) +
						 u_old(x, y - 1) - d * d * f(x, y)) /
					    4.0;
				}
			}

			double max_residual = 0.0;
			for (int y = 1; y < n - 1; y++) {
				for (int x = 1; x < n - 1; x++) {
					double laplacian = five_point_laplacian(
					    u_new(x, y),
					    u_new(x - 1, y),
					    u_new(x + 1, y),
					    u_new(x, y - 1),
					    u_new(x, y + 1),
					    d
					);
					double residual = std::abs(f(x, y) - laplacian);
					if (residual > max_residual)
						max_residual = residual;
				}
			}
			if (max_residual < residual_tolerance) {
				std::cout << "Residual " << max_residual << " reached after "
				          << i + 1 << " iterations\n";
				break;
			}
		}

		double max_error = 0.0;
		for (int y = 1; y < n - 1; y++) {
			for (int x = 1; x < n - 1; x++) {
				double error = std::abs(
				    u_new(x, y) - exact_solution(min + x * d, min + y * d)
				);
				if (error > max_error)
					max_error = error;
			}
		}

		std::cout << "Max error: " << max_error << '\n';
		auto end = std::chrono::steady_clock::now();
		double elapsed = std::chrono::duration<double>(end - start).count();
		std::cout << "Time: " << elapsed << " s\n";
	}
}

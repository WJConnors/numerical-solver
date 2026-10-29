#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <utility>

#include "finite_difference.hpp"
#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

// Grid hierarchy: 3, 5, 9, 17, 33, 65, 129 etc..
constexpr int finest = 129;

constexpr double min = 0.0;
constexpr double max = 1.0;

constexpr double weighting = 2.0 / 3.0;

constexpr double residual_tolerance = 1e-8;

// To be used for final verification, assuming u is unknown otherwise
// u(x,y) = sin(pi.x).sin(pi.y)
double exact_solution(double x, double y) { return std::sin(pi * x) * std::sin(pi * y); }
// Laplacian
//  f(x,y) = -2.pi^2.sin(pi.x).sin(pi.y)
double exact_laplacian(double x, double y) { return -2 * pi * pi * exact_solution(x, y); }

bool is_even(int i) { return i % 2 == 0; }

void restrict(const Grid2D &fine_grid, Grid2D &coarse_grid) {
	for (int y = 1; y < coarse_grid.ny - 1; y++) {
		for (int x = 1; x < coarse_grid.nx - 1; x++) {
			int centre_x = x * 2;
			int centre_y = y * 2;

			coarse_grid(x, y) =
			    (1 * fine_grid(centre_x - 1, centre_y - 1) + 2 * fine_grid(centre_x, centre_y - 1) +
				 1 * fine_grid(centre_x + 1, centre_y - 1) + 2 * fine_grid(centre_x - 1, centre_y) +
				 4 * fine_grid(centre_x, centre_y) + 2 * fine_grid(centre_x + 1, centre_y) +
				 1 * fine_grid(centre_x - 1, centre_y + 1) + 2 * fine_grid(centre_x, centre_y + 1) +
				 1 * fine_grid(centre_x + 1, centre_y + 1)) /
			    16.0;
		}
	}
}

void prolong(const Grid2D &coarse_grid, Grid2D &fine_grid) {
	for (int y = 1; y < fine_grid.ny - 1; y++) {
		for (int x = 1; x < fine_grid.nx - 1; x++) {
			if (is_even(x) && is_even(y)) {
				fine_grid(x, y) = coarse_grid(x / 2, y / 2);
			} else if (!is_even(x) && is_even(y)) {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) = (coarse_grid(cx, cy) + coarse_grid(cx + 1, cy)) / 2.0;
			} else if (is_even(x) && !is_even(y)) {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) = (coarse_grid(cx, cy) + coarse_grid(cx, cy + 1)) / 2.0;
			} else {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) = (coarse_grid(cx, cy) + coarse_grid(cx + 1, cy) +
				                   coarse_grid(cx, cy + 1) + coarse_grid(cx + 1, cy + 1)) /
				                  4.0;
			}
		}
	}
}

void weighted_jacobi(Grid2D &u, const Grid2D &f, double d, int sweeps) {
	Grid2D u_old = make_grid2d(u.nx, u.ny);
	for (int i = 0; i < sweeps; i++) {
		std::swap(u_old, u);
		for (int y = 1; y < u_old.ny - 1; y++) {
			for (int x = 1; x < u_old.nx - 1; x++) {
				double jacobi_value = (u_old(x + 1, y) + u_old(x - 1, y) + u_old(x, y + 1) +
				                       u_old(x, y - 1) - d * d * f(x, y)) /
				                      4.0;
				u(x, y) = (1 - weighting) * u_old(x, y) + weighting * jacobi_value;
			}
		}
	}
}

void calculate_residual(const Grid2D &u, const Grid2D &f, Grid2D &residual, double d) {
	for (int y = 1; y < residual.ny - 1; y++) {
		for (int x = 1; x < residual.nx - 1; x++) {
			double laplacian = five_point_laplacian(
			    u(x, y),
			    u(x - 1, y),
			    u(x + 1, y),
			    u(x, y - 1),
			    u(x, y + 1),
			    d
			);
			residual(x, y) = f(x, y) - laplacian;
		}
	}
}

double abs_max_grid_value(const Grid2D &grid) {
	double max_value = 0.0;
	for (int y = 0; y < grid.ny; y++) {
		for (int x = 0; x < grid.nx; x++) {
			double value = std::abs(grid(x, y));
			if (value > max_value)
				max_value = value;
		}
	}
	return max_value;
}

void v_cycle(Grid2D &u, const Grid2D &f, double d) {
	int n = u.nx;

	if (n == 3) {
		u(1, 1) = -d * d * f(1, 1) / 4.0;
		return;
	}

	int coarse_n = (n + 1) / 2;

	Grid2D fine_residual = make_grid2d(n, n);
	Grid2D coarse_residual = make_grid2d(coarse_n, coarse_n);
	Grid2D fine_correction = make_grid2d(n, n);
	Grid2D coarse_correction = make_grid2d(coarse_n, coarse_n);

	// Pre-smoothing
	weighted_jacobi(u, f, d, 3);
	// Current residual
	calculate_residual(u, f, fine_residual, d);
	restrict(fine_residual, coarse_residual);

	v_cycle(coarse_correction, coarse_residual, 2.0 * d);
	// Correct for error
	prolong(coarse_correction, fine_correction);
	for (int y = 1; y < n - 1; y++) {
		for (int x = 1; x < n - 1; x++) {
			u(x, y) += fine_correction(x, y);
		}
	}

	// Post-smoothing
	weighted_jacobi(u, f, d, 3);
}

int main() {
	double d = (max - min) / (finest - 1);
	Grid2D f = make_grid2d(finest, finest);
	for (int y = 0; y < finest; y++) {
		for (int x = 0; x < finest; x++) {
			f(x, y) = exact_laplacian(min + x * d, min + y * d);
		}
	}

	Grid2D u = make_grid2d(finest, finest);
	Grid2D residual = make_grid2d(finest, finest);
	calculate_residual(u, f, residual, d);
	double initial_residual = abs_max_grid_value(residual);
	double current_residual = initial_residual;

	auto start = std::chrono::steady_clock::now();
	while (current_residual / initial_residual > residual_tolerance) {
		v_cycle(u, f, d);
		calculate_residual(u, f, residual, d);
		current_residual = abs_max_grid_value(residual);
	}
	auto end = std::chrono::steady_clock::now();
	double elapsed = std::chrono::duration<double>(end - start).count();

	double max_error = 0.0;
	for (int y = 0; y < u.ny; y++) {
		for (int x = 0; x < u.nx; x++) {
			double exact = exact_solution(min + x * d, min + y * d);
			double error = std::abs(u(x, y) - exact);

			if (error > max_error)
				max_error = error;
		}
	}

	std::cout << "Solved in " << elapsed << "s\nMax residual: " << current_residual
	          << "\nMax error: " << max_error << '\n';
}

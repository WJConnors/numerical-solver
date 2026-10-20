#include <cmath>
#include <iostream>
#include <numbers>
#include <utility>

#include "finite_difference.hpp"
#include "grid2d.hpp"

constexpr double pi = std::numbers::pi;

constexpr int fine = 65;
constexpr int coarse = 33;

constexpr double min = 0.0;
constexpr double max = 1.0;

constexpr double weighting = 2.0 / 3.0;

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

bool is_even(int i) { return i % 2 == 0; }

void restrict(const Grid2D &fine_grid, Grid2D &coarse_grid) {
	for (int y = 1; y < coarse_grid.ny - 1; y++) {
		for (int x = 1; x < coarse_grid.nx - 1; x++) {
			int centre_x = x * 2;
			int centre_y = y * 2;

			coarse_grid(x, y) = (1 * fine_grid(centre_x - 1, centre_y - 1) +
			                     2 * fine_grid(centre_x, centre_y - 1) +
			                     1 * fine_grid(centre_x + 1, centre_y - 1) +
			                     2 * fine_grid(centre_x - 1, centre_y) +
			                     4 * fine_grid(centre_x, centre_y) +
			                     2 * fine_grid(centre_x + 1, centre_y) +
			                     1 * fine_grid(centre_x - 1, centre_y + 1) +
			                     2 * fine_grid(centre_x, centre_y + 1) +
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
				fine_grid(x, y) =
				    (coarse_grid(cx, cy) + coarse_grid(cx + 1, cy)) / 2.0;
			} else if (is_even(x) && !is_even(y)) {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) =
				    (coarse_grid(cx, cy) + coarse_grid(cx, cy + 1)) / 2.0;
			} else {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) =
				    (coarse_grid(cx, cy) + coarse_grid(cx + 1, cy) +
					 coarse_grid(cx, cy + 1) + coarse_grid(cx + 1, cy + 1)) /
				    4.0;
			}
		}
	}
}

void weighted_jacobi(Grid2D &u_old, Grid2D &u_new, const Grid2D &f, double d) {
	std::swap(u_old, u_new);
	for (int y = 1; y < u_old.ny - 1; y++) {
		for (int x = 1; x < u_old.nx - 1; x++) {
			double jacobi_value =
			    (u_old(x + 1, y) + u_old(x - 1, y) + u_old(x, y + 1) +
				 u_old(x, y - 1) - d * d * f(x, y)) /
			    4.0;
			u_new(x, y) =
			    (1 - weighting) * u_old(x, y) + weighting * jacobi_value;
		}
	}
}

void calculate_residual(
    const Grid2D &u,
    const Grid2D &f,
    Grid2D &residual,
    double d
) {
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

int main() {
	double fine_d = (max - min) / (fine - 1);
	double coarse_d = (max - min) / (coarse - 1);
	Grid2D f = make_grid2d(fine, fine);
	for (int y = 0; y < fine; y++) {
		for (int x = 0; x < fine; x++) {
			f(x, y) = exact_laplacian(min + x * fine_d, min + y * fine_d);
		}
	}

	Grid2D u_old = make_grid2d(fine, fine);
	Grid2D u_new = make_grid2d(fine, fine);
	Grid2D fine_residual = make_grid2d(fine, fine);
	Grid2D coarse_residual = make_grid2d(coarse, coarse);
	Grid2D coarse_correction_old = make_grid2d(coarse, coarse);
	Grid2D coarse_correction_new = make_grid2d(coarse, coarse);
	Grid2D fine_correction = make_grid2d(fine, fine);

	for (int i = 0; i < 3; i++) {
		weighted_jacobi(u_old, u_new, f, fine_d);
	}

	calculate_residual(u_new, f, fine_residual, fine_d);
	restrict(fine_residual, coarse_residual);

	for (int i = 0; i < 1000; i++) {
		weighted_jacobi(
		    coarse_correction_old,
		    coarse_correction_new,
		    coarse_residual,
		    coarse_d
		);
	}

	prolong(coarse_correction_new, fine_correction);
}

#include <iostream>

#include "grid2d.hpp"

constexpr int fine = 65;
constexpr int coarse = 33;

bool is_even(int i) {
	return i % 2 == 0;
}

void restrict(const Grid2D& fine_grid, Grid2D& coarse_grid) {
	for (int y = 1; y < coarse_grid.ny - 1; y++) {
		for (int x = 1; x < coarse_grid.nx - 1; x++) {
			int centre_x = x * 2;
			int centre_y = y * 2;

			coarse_grid(x, y) = (
					1 * fine_grid(centre_x - 1, centre_y - 1) +
					2 * fine_grid(centre_x, centre_y - 1) +
					1 * fine_grid(centre_x + 1, centre_y - 1) +
					2 * fine_grid(centre_x - 1, centre_y) + 
					4 * fine_grid(centre_x, centre_y) +
					2 * fine_grid(centre_x + 1, centre_y) +
					1 * fine_grid(centre_x - 1, centre_y + 1) +
					2 * fine_grid(centre_x, centre_y + 1) +
					1 * fine_grid(centre_x + 1, centre_y + 1)
				) / 16.0;
		}
	}
}

void prolong(const Grid2D& coarse_grid, Grid2D& fine_grid) {
	for (int y = 1; y < fine_grid.ny - 1; y++) {
		for (int x = 1; x < fine_grid.nx - 1; x++){
			if (is_even(x) && is_even(y)) {
				fine_grid(x, y) = coarse_grid(x / 2, y / 2);
				}
			else if (!is_even(x) && is_even(y)) {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) = 
					(coarse_grid(cx, cy)
					+ coarse_grid(cx + 1, cy))
					/ 2.0;
			}
			else if (is_even(x) && !is_even(y)){
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) =
					(coarse_grid(cx, cy)
					+ coarse_grid(cx, cy + 1))
					/ 2.0;
			}
			else {
				int cx = x / 2;
				int cy = y / 2;
				fine_grid(x, y) = 
				(coarse_grid(cx, cy)
				+ coarse_grid(cx + 1, cy)
				+ coarse_grid(cx, cy + 1)
				+ coarse_grid(cx + 1, cy + 1))
				/ 4.0;
			}
		}
	}
}

int main() {
	std::cout << "Hello, multigrid\n";
}

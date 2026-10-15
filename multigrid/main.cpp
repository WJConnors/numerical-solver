#include <iostream>

#include "grid2d.hpp"

constexpr int fine = 65;
constexpr int coarse = 33;

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

int main() {
	std::cout << "Hello, multigrid\n";
}

#pragma once

#include <vector>

struct Grid2D {
	int nx;
	int ny;
	std::vector<double> values;

	double& operator()(int x, int y) {
		return values[y * nx + x];
	}	
};

inline Grid2D make_grid2d(int nx, int ny) {
	Grid2D grid;

	grid.nx = nx;
	grid.ny = ny;
	grid.values = std::vector<double>(nx * ny);

	return grid;
}

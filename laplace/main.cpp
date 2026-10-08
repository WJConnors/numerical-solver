#include <iostream>
#include <numbers>
#include <cmath>

#include "grid2d.hpp"
#include "finite_difference.hpp"

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

double five_point_laplacian(
	double centre,
	double x_prev,
	double x_next,
	double y_prev,
	double y_next,
	double d
) {
	return finite_difference::second_derivative(x_prev, centre, x_next, d)
		+ finite_difference::second_derivative(y_prev, centre, y_next, d);
}

int main() {

	for (int n : ns) {
		std::cout << "n = " << n << '\n';
		
		double d = (max - min) / (n -1);

		Grid2D values = make_grid2d(n, n);

		for (int y = 0; y < n; y++) {
			for (int x = 0; x < n; x++) {
				values(x, y) = exact_solution(
					min + x * d,
					min + y * d
				);
				
			}
		}

		Grid2D fp_laplacian_grid = make_grid2d(n, n);

		double max_error = 0.0;

		for (int y = 1; y < n - 1; y++) {
			for (int x = 1; x < n - 1; x++) {
				double exact = exact_laplacian(
					min + x * d,
					min + y * d
				);
				
				fp_laplacian_grid(x, y) = five_point_laplacian(
					values(x, y),
					values(x - 1, y),
					values(x + 1, y),
					values(x, y - 1),
					values(x, y + 1),
					d
				);

				double error = std::abs(fp_laplacian_grid(x, y) - exact);
				if (error > max_error) max_error = error;
			}
		}

		std::cout << "Max Error = " << max_error << "\n\n";
	}
}

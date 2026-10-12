#pragma once

inline double finite_difference_first_derivative(double prev, double next, double d) {
	return (next - prev) / (2.0 * d);
}

inline double finite_difference_second_derivative(double prev, double cur, double next, double d) {
	return (next - 2.0 * cur + prev) / (d * d);
}

inline double five_point_laplacian(
	double centre,
	double x_prev,
	double x_next,
	double y_prev,
	double y_next,
	double d
) {
	return finite_difference_second_derivative(x_prev, centre, x_next, d)
		+ finite_difference_second_derivative(y_prev, centre, y_next, d);
}

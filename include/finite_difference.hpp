#pragma once

inline double finite_difference_first_derivative(double prev, double next, double d) {
	return (next - prev) / (2.0 * d);
}

inline double finite_difference_second_derivative(double prev, double cur, double next, double d) {
	return (next - 2.0 * cur + prev) / (d * d);
}

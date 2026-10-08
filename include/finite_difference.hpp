#pragma once

namespace finite_difference {

inline double first_derivative(double prev, double next, double d) {
	return (next - prev) / (2.0 * d);
}

inline double second_derivative(double prev, double cur, double next, double d) {
	return (next - 2.0 * cur + prev) / (d * d);
}

}

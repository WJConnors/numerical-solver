#include <iostream>
#include <cmath>
#include <numbers>
#include <vector>

constexpr double pi = std::numbers::pi;

//Interval: 0 <= x <= 1
constexpr int n = 11;
constexpr double lower_limit = 0.0;
constexpr double upper_limit = 1.0;

double u(double x) {
	return std::sin(x * pi);
}

double exact_first_derivative(double x) {
	return pi * std::cos(x * pi);
}

double exact_second_derivative(double x) {
	return -pi * pi
		* std::sin(x * pi);
}

int main() {
    /*std::cout << "u(0) = " << u(0) << std::endl;
    std::cout << "u(0.5) = " << u(0.5) << std::endl;
    std::cout << "u(1) = " << u(1) << std::endl;

    std::cout << "u'(0) = " << exact_first_derivative(0) << std::endl;
    std::cout << "u'(0.5) = " << exact_first_derivative(0.5) << std::endl;
    std::cout << "u'(1) = " << exact_first_derivative(1) << std::endl;

    std::cout << "u''(0) = " << exact_second_derivative(0) << std::endl;
    std::cout << "u''(0.5) = " << exact_second_derivative(0.5) << std::endl;
    std::cout << "u''(1) = " << exact_second_derivative(1) << std::endl;*/

    std::vector<double> x(n), values(n);
    double dx = (upper_limit - lower_limit) / (n - 1);
    for (int i = 0; i < n; i++) {
    	x[i] = lower_limit + i * dx;
    	std::cout << "x = " << x[i] << std::endl;
    	values[i] = u(x[i]);
    	std::cout << "u(x) = " << values[i] << std::endl;
    }
}

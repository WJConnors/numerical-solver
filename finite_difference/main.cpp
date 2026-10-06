#include <iostream>
#include <cmath>
#include <numbers>

constexpr double pi = std::numbers::pi;

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
    std::cout << "u(0) = " << u(0) << std::endl;
    std::cout << "u(0.5) = " << u(0.5) << std::endl;
    std::cout << "u(1) = " << u(1) << std::endl;

    std::cout << "u'(0) = " << exact_first_derivative(0) << std::endl;
    std::cout << "u'(0.5) = " << exact_first_derivative(0.5) << std::endl;
    std::cout << "u'(1) = " << exact_first_derivative(1) << std::endl;

    std::cout << "u''(0) = " << exact_second_derivative(0) << std::endl;
    std::cout << "u''(0.5) = " << exact_second_derivative(0.5) << std::endl;
    std::cout << "u''(1) = " << exact_second_derivative(1) << std::endl;
}

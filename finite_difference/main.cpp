#include <iostream>
#include <cmath>
#include <numbers>
#include <vector>

constexpr double pi = std::numbers::pi;

//Interval: 0 <= x <= 1
//Doubling n reduces max error by ~75%
constexpr int ns[] = {11, 21, 41, 81};
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
	for (int n : ns) {

		std::vector<double> x(n), values(n);
		double dx = (upper_limit - lower_limit) / (n - 1);

		std::cout << '\n' << "n = " << n << " dx = " << dx << '\n';


	    //std::cout << "Original Values" << '\n';
	    for (int i = 0; i < n; i++) {
	    	x[i] = lower_limit + i * dx;
	    	//std::cout << "x = " << x[i] << '\n';
	    	values[i] = u(x[i]);
	    	//std::cout << "u(x) = " << values[i] << '\n';
	    }
	    
	    //std::cout << std::endl;
	    std::cout << "Calculating First Derivative" << '\n';

	    std::vector<double> d1(n), exact_d1(n);
	    double max_error = 0.0;
	    for (int i = 1; i < n - 1; i++) {
	    	exact_d1[i] = exact_first_derivative(x[i]);
	    	d1[i] = (values[i+1] - values[i-1])
	    		/ (2.0 * dx);
	    	double error = std::abs(d1[i]- exact_d1[i]);
	    	if (error > max_error) max_error = error;
	    	/*std::cout << "x = " << x[i]
	    		<< " u'(x) = " << d1[i]
	    		<< " exact u'(x) " << exact_d1[i]
	    		<< " error = " << error
	    		<< '\n';*/
	    }
	    std::cout << "Max Error = " << max_error << '\n';

	    std::cout << '\n';
		std::cout << "Calculating Second Derivative" << '\n';
	    
	    std::vector<double> d2(n), exact_d2(n);
	    max_error = 0.0;
	    for (int i = 1; i < n - 1; i++) {
	    	exact_d2[i] = exact_second_derivative(x[i]);
	    	d2[i] = (values[i+1] - 2.0 * values[i] + values[i-1])
	    		/ (dx * dx);
	    	double error = std::abs(d2[i]- exact_d2[i]);
	    	if (error > max_error) max_error = error;
	    	/*std::cout << "x = " << x[i]
	    	    	<< " u''(x) = " << d2[i]
	    	    	<< " exact u''(x) " << exact_d2[i]
	    	    	<< " error = " << error
	    	    	<< '\n';*/
	    }
	    std::cout << "Max Error = " << max_error << '\n';
    }
    
}

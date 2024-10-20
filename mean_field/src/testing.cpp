#include <iostream>
/* #include <omp.h> */
#include "DESolver.h"
#include <unordered_map>

int main (int argc, char *argv[]) {
	DESolver desolver(0);
	std::vector<double> y = {40.1, 1.1};
	std::unordered_map<std::string, double> p;
	p["us"] = 40;
	p["k0"] = 0.0125;
	p["k1"] = 0.0111;
	p["alpha"] = 1;
	p["n"] = 0.5;
	desolver.set_params(p);
	if (!desolver.check_param_region()) {
		std::cerr << "Error: Parameters are not in the correct region" << std::endl;
		return 1;
	}
	desolver.set_y(y[0], y[1]);
	desolver.initialize();
	double period = desolver.get_period();
	std::cout << "Period: " << period << std::endl;
	double t = 0;
	while (t < 500*period) {
		desolver.solve(t, t+period/100);
		t = t+period/100;
	}
	y = desolver.get_y();
	std::cout << "Density values: ";
	std::cout << t << " " << y[0] << " " << y[1] << std::endl;
	return 0;
}

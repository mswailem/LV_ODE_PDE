#include <iostream>
/* #include <omp.h> */
#include "DESolver.h"
#include <unordered_map>

int main (int argc, char *argv[]) {

	DESolver desolver(0);
	std::vector<double> y = {0.5, 0.5};
	std::unordered_map<std::string, double> p;
	p["us"] = 1;
	p["k0"] = 0.45;
	p["k1"] = 0.1;
	p["alpha"] = 1;
	p["n"] = 1;
	double t0 = 0;
	double tf = 10;
	double dt = 1;
	desolver.set_params(p);
	if (!desolver.check_param_region()) {
		std::cerr << "Error: Parameters are not in the correct region" << std::endl;
		return 1;
	}
	desolver.set_y(y[0], y[1]);
	desolver.initialize();
	std::cout << "Parameters: ";
	desolver.d_print_params();
	std::cout << "Omega0: ";
	desolver.d_print_omega0();
	std::cout << "Lambda0: ";
	desolver.d_print_lambda0();
	std::cout << "Period: ";
	desolver.d_print_period();

	for (double t = t0; t < tf; t += dt) {
		std::cout << "Density values: ";
		std::cout << t << " " << y[0] << " " << y[1] << std::endl;
		std::cout << "K: ";
		desolver.d_print_k(t);
		std::cout << "Lambda: ";
		desolver.d_print_lambda(t);
		desolver.solve(t, t+dt);
		y = desolver.get_y();
	}
	std::cout << "Density values: ";
	std::cout << tf << " " << y[0] << " " << y[1] << std::endl;
	std::cout << "K: ";
	desolver.d_print_k(tf);
	std::cout << "Lambda: ";
	desolver.d_print_lambda(tf);
	return 0;
}

#include <iostream>
#include <vector>
#include "DESolver.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

int main (int argc, char *argv[]) {
	DESolver solver(1, 1, 0.22, 0, 1, 1, 1.1, 1.1, 0);
	solver.initialize();
	double dt = 0.01;
	std::cout << solver.get_y()[0] << " " << solver.get_y()[1] << std::endl;
	solver.solve(0, 1, dt);
	std::cout << solver.get_y()[0] << " " << solver.get_y()[1] << std::endl;
	solver.solve(0, 1, dt);
	std::cout << solver.get_y()[0] << " " << solver.get_y()[1] << std::endl;
	DESolver solver2(1, 1, 0.22, 0, 1, 1, 1.1, 1.1, 0);
	std::cout << solver2.get_y()[0] << " " << solver2.get_y()[1] << std::endl;
	solver2.solve(0, 2, dt);
	std::cout << solver2.get_y()[0] << " " << solver2.get_y()[1] << std::endl;
	return 0;
}

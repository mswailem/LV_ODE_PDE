#include <iostream>
#include <omp.h>
#include "phase.h"
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "DESolver.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


// NOTE: Currently testing the phase plot of k1 vs alpha, will end up refactoring this better into the code 

int main (int argc, char *argv[]) {

	int t0 = std::stoi(argv[1]);
	int points = std::stoi(argv[2]);
	double ustar = std::stod(argv[3]);
	double vstar = std::stod(argv[4]);
	double k0 = std::stod(argv[5]);
	double n = std::stod(argv[6]);
	double step = std::stod(argv[7]);

	double alpha = 1;
	double k1;

	// For parallel processing of k1 loop, we need to iterate over integers
	int num_steps = static_cast<int>((k0 / step) + 1);

	std::string command = "mkdir -p ../output/testing";
	std::system(command.c_str());
	std::ofstream out_file;
	out_file.open("../output/testing/phase_k1_vs_alpha.dat");

	std::vector<std::pair<double, double>> fps;
	fps.clear(); // Just in case
	
	int num_of_points;

	while (alpha >= 0) {
		#pragma omp parallel for private(k1) shared(out_file, alpha)
		for (int i = 0; i < num_steps; i++) {

			k1 = i * step;
			DESolver desolver(0);

			desolver.set_k1(k1);
			desolver.set_n(n);
			desolver.set_ustar(ustar);
			desolver.set_vstar(vstar);
			desolver.set_k0(k0);
			desolver.set_alpha(alpha);
			desolver.set_y(ustar*(1+pow(10,-3)), vstar*(1+pow(10,-3)));

			fps = compute_fixed_points(desolver, t0, points, 1e-4);
			num_of_points = fps.size();

			std::cout << alpha << " " << k1 << std::endl;
			#pragma omp critical // Protects file writing 
			{
			out_file << alpha << " " << k1 << " " << num_of_points << std::endl;
			}

			k1 += step;
		}
		alpha -= step;
	}

	return 0;
}

#include <iostream>
#include <omp.h>
#include "programs.h"
#include "theory.h"
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "DESolver.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


// TODO: Will need to move this function somewhere else after refactoring the code

void phase_space(int argc, char *argv[]){
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
			desolver.set_fp(ustar, vstar);
			desolver.set_k0(k0);
			desolver.set_alpha(alpha);
			desolver.set_y(ustar*(1+pow(10,-3)), vstar*(1+pow(10,-3)));

			fps = compute_stationary_points(desolver, t0, points, 1e-4);
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
}

// TODO: Will need to move this function somewhere else after refactoring the code

void dispersion_relation(int argc, char *argv[]) {
	double du = std::stod(argv[1]);
	double dv = std::stod(argv[2]);

	double u0 = 1;
	double v0 = 1;
	double k0 = 0.245;
	double wn = 0;

	int points = 20;
	double wnf = 2;
	double step = (wnf-wn)/points;

	std::string command = "mkdir -p ../output/testing/eigenvalues";
	std::system(command.c_str());
	std::ofstream out_file;
	out_file.open("../output/testing/eigenvalues/" + std::string(argv[1]) + "_" + std::string(argv[2]) + ".dat");

	#pragma omp parallel for private(wn) shared(out_file)
	for (int i = 0; i < points; i++) {

		wn = i * step;
		gsl_vector_complex *eigenvalues = get_eigenvalues(du, dv, u0, v0, k0, wn);

		gsl_complex lambda1 = gsl_vector_complex_get(eigenvalues, 0);
		gsl_complex lambda2 = gsl_vector_complex_get(eigenvalues, 1);

		// TODO: change writing data as csv files so that it is easier to handle and more transparent
		#pragma omp critical
		{
			out_file << wn << " " << GSL_REAL(lambda1) << " " << GSL_IMAG(lambda1) << " " << GSL_REAL(lambda2) << " " << GSL_IMAG(lambda2) << std::endl;
		}

		gsl_vector_complex_free(eigenvalues);
	}
}

// NOTE: I am currently implementing the function to get the dispersion relation for the system, will need to move this after refactoring

int main (int argc, char *argv[]) {
	
	std::unordered_map<std::string, double> p;
	p["us"] = 1;
	p["vs"] = 1;
	p["k0"] = 0.245;
	p["k1"] = 0;
	p["n"] = 0;
	p["alpha"] = 1;
	VaryingParam k1 = VaryingParam("k1", 0, 0.245, 0.001);
	VaryingParam n = VaryingParam("n", 0, 0.57, 0.001);
	stability(p, n, k1, 1000);

	return 0;
}

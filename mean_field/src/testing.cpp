#include <iostream>
#include <omp.h>
#include "iohandler.h"
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
// NOTE: I am leaving this here for the parallel processing so that I can revisit it and implement it later in my programs.h file
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

// NOTE: I am currently implementing the function to get the dispersion relation for the system, will need to move this after refactoring

int main (int argc, char *argv[]) {
	
	std::string picked_program_name = fzf_pick(get_program_names(), "Pick program: ");
	std::cout << "Picked program: " << picked_program_name << std::endl;
	Program picked_program = get_program(picked_program_name);
	std::vector<VaryingParam> varying_params = pick_varying_params(picked_program);
	for (auto const& x : varying_params) {
		std::cout << x.name << " " << x.start << " " << x.end << " " << x.step << std::endl;
	}
	return 0;
}

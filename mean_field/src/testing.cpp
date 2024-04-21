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

// NOTE: Might change this so that instead each program has a function member variable
void run_program(Program program, std::unordered_map<std::string, double> params, std::vector<VaryingParam> vs, std::string filename) {
	if (program.name == "bifurcation diagram") {
		if (vs.size() == 1) {
			bifurcation_diagram(params, vs[0], params["t0"], params["points_in_period"], filename);
		} else if (vs.size() == 2) {
			bifurcation_diagram(params, vs[0], vs[1], params["t0"], params["points_in_period"], filename);
		} else {
			throw std::invalid_argument("Invalid number of varying parameters");
		}
	} else if (program.name == "time series") {
		time_series(params, params["t0"], params["tf"], params["dt"], {params["a0"], params["b0"]}, filename);
	} else if (program.name == "stability") {
		stability(params, vs[0], vs[1], params["points_in_period"], filename);
	} else if (program.name == "dispersion relation") {
		dispersion_relation(params, vs[0], filename);
	} else {
		throw std::invalid_argument("Invalid program name");
	}
}

// TODO: Move this to the main.cpp file, and remove the unnecessary includes from all the files
int main (int argc, char *argv[]) {
	
	std::string picked_program_name = fzf_pick(get_program_names(), "Pick program: ");
	Program picked_program = get_program(picked_program_name);
	std::vector<VaryingParam> varying_params = pick_varying_params(picked_program);
	std::unordered_map<std::string, double> fixed_params = get_params(picked_program, varying_params);
	std::string filename = get_filename();
	run_program(picked_program, fixed_params, varying_params, filename);
	return 0;
}

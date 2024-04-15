#include <iostream>
#include <vector>
#include "DESolver.h"
#include <string>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <cstdlib>
#include <cmath>
#include "programs.h"

// TODO: Implement the bifurcation diagram for different vairables, and make the pipline more efficient
// TODO: Maybe move some of these functions to a separate file?
// TODO: Create a function that automatically checks the inputs, and maybe merge this with printUsage?
// TODO: The parameters of the system can maybe be handeled by a struct instead?

double t0, tf, dt, ustar, vstar, k0, k1, n, alpha, step, wavenumber;
int t_points, points, type;
std::vector<double> y = {0.0, 0.0};
std::string output_type;

//Prints the usage of the program
void printUsage(bool failed) { 
    if (output_type == "time_series" && failed) {
		std::cerr << "Required arguments for time_series: t0 tf dt ustar vstar k0 k1 n alpha initial_conditions" << std::endl;
	} else if (output_type == "stability_k1_vs_n" && failed) {
		std::cerr << "Required arguments for stability diagram as a function of k1 and n: points k0 step ustar vstar wavenumber" << std::endl;
	} else if (output_type == "stability_ustar_vs_vstar" && failed) {
		std::cerr << "Required arguments for stability diagram as a function of ustar and vstar: points n step" << std::endl;
	} else if (output_type == "bifurcation_alpha" && failed) {
		std::cerr << "Required arguments for bifurcation diagram as a function of alpha: t0 points ustar vstar k0 k1 n step" << std::endl;
	} else {
		std::cerr << "Program type not found, available options are: time_series, stability_k1_vs_n, stability_ustar_vs_vstar, bifurcation_alpha" << std::endl;
	}
}

int main(int argc, char* argv[]) {
    if (argc > 1 && argv[1] != nullptr) {
		output_type = argv[1];
		if (output_type == "time_series") {
			/* run_time_series(argc, argv); */
		} else if (output_type == "stability_k1_vs_n") {
			/* run_stability(argc,argv, "k1_vs_n"); */
		} else if (output_type == "stability_ustar_vs_vstar") {
			/* run_stability(argc,argv, "ustar_vs_vstar"); */
		} else if (output_type == "bifurcation_alpha") {
			/* run_bifurcation(argc, argv, "alpha"); */
			return 0;
		}
    } else {
        printUsage(true);
    }
}

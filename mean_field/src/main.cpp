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

double t0;
double tf;
double dt;
int t_points;
double ustar;
double vstar;
double k0;
double k1;
double n;
double alpha;
int points;
double step;
std::vector<double> y = {0.0, 0.0};
int type;
std::string output_type;

void printUsage(bool failed) { //Prints the usage of the program
    if (output_type == "time_series" && failed) {
		std::cerr << "Required arguments for time_series: t0 tf dt ustar vstar k0 k1 n alpha initial_conditions" << std::endl;
	} else if (output_type == "stability_k1_vs_n" && failed) {
		std::cerr << "Required arguments for stability diagram as a function of k1 and n: points k0 step ustar vstar" << std::endl;
	} else if (output_type == "stability_ustar_vs_vstar" && failed) {
		std::cerr << "Required arguments for stability diagram as a function of ustar and vstar: points n step" << std::endl;
	} else if (output_type == "bifurcation_alpha" && failed) {
		std::cerr << "Required arguments for bifurcation diagram as a function of alpha: t0 points ustar vstar k0 k1 n step" << std::endl;
	} else {
		std::cerr << "Program type not found, available options are: time_series, stability_k1_vs_n, stability_ustar_vs_vstar, bifurcation_alpha" << std::endl;
	}
}

gsl_matrix* compute_fm(DESolver& solver) { //This computes the fundemental matrix for the system
		// Variable defintions
		tf = solver.get_period();
		dt = tf/points;

		//Allocating memory for gsl variables
		static gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);

		//Solving the system for the the (1,0) vector
		solver.set_y(1, 0);
		solver.solve(0, tf, dt);
		y = solver.get_y();

		//Assigning the values of the fundemental matrix
		gsl_matrix_set(fundemental_matrix, 0, 0, y[0]);
		gsl_matrix_set(fundemental_matrix, 0, 1, y[1]);

		//Solving the system for the (0,1) vector
		solver.set_y(0, 1);
		solver.solve(0, tf, dt);
		y = solver.get_y();

		//Assigning the values of the fundemental matrix
		gsl_matrix_set(fundemental_matrix, 1, 0, y[0]);
		gsl_matrix_set(fundemental_matrix, 1, 1, y[1]);


		return fundemental_matrix;
		gsl_matrix_free(fundemental_matrix);
}

void run_time_series(int argc, char* argv[]) {

	//Check arguments are correct
	if (argc < 13) {
		printUsage(true);
		return;
	}
	
	//Parse arguments
	t0 = std::stod(argv[2]);
	tf = std::stod(argv[3]);
	dt = std::stod(argv[4]);
	ustar = std::stod(argv[5]);
	vstar = std::stod(argv[6]);
	k0 = std::stod(argv[7]);
	k1 = std::stod(argv[8]);
	n = std::stod(argv[9]);
	alpha = std::stod(argv[10]);
	y[0] = std::stod(argv[11]);
	y[1] = std::stod(argv[12]);

	//Start the simulation
	DESolver desolver(ustar, vstar, k0, k1, n, alpha, y[0], y[1], 0);

	//Output file named by te parameters
	std::ofstream out_file("../output/time_series/" + std::string(argv[5]) + "_" + std::string(argv[6]) + "_" + std::string(argv[7]) + "_" + std::string(argv[8]) + "_" + std::string(argv[9]) + "_" + std::string(argv[10]) + ".dat");
	//Loop over time and solve the system at each time step while printing output
	for (double t = t0; t < tf; t += dt) {
		out_file << t << " " << y[0] << " " << y[1] << std::endl;
		desolver.solve(t, t+dt, dt);
		y = desolver.get_y();
	}
	out_file << tf << " " << y[0] << " " << y[1] << std::endl;
}

void run_stability(int argc, char* argv[], std::string vars) {

	//Check arguments are correct (might move this to its own function
	if (vars == "k1_vs_n") {
		if (argc < 7) {
			printUsage(true);
			return;
		} else {
			//Parse arguments
			points = std::stoi(argv[2]);
			k0 = std::stod(argv[3]);
			step = std::stod(argv[4]);
			ustar = std::stod(argv[5]);
			vstar = std::stod(argv[6]);
			n = step;
		}
	} else if (vars == "ustar_vs_vstar") {
		if (argc < 5) {
			printUsage(true);
			return;
		} else {
			//Parse arguments
			points = std::stoi(argv[2]);
			n = std::stod(argv[3]);
			step = std::stod(argv[4]);
			ustar = step;
		}
	}
	
	//Allocating memory for gsl variables
	gsl_vector_complex *floquet_multipliers = gsl_vector_complex_alloc(2);
	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
	gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);

	//Defining the linear DESolver (the parameters will be specificed later)
	DESolver desolver(1);

	std::string command = "mkdir -p ../output/fm/" + vars;
	std::system(command.c_str());
	std::ofstream out_file;
	if (vars == "k1_vs_n") {
		out_file.open("../output/fm/" + vars + "/k0=" + std::string(argv[3]) + "_ustar=" + std::string(argv[5]) + "_vtar=" + std::string(argv[6]) + ".dat");
	} else if (vars == "ustar_vs_vstar") {
		out_file.open("../output/fm/" + vars + "/n=" + std::string(argv[3]) + ".dat");
	}

	while ((n < 2.1 && vars == "k1_vs_n") || (ustar <= 5 && vars == "ustar_vs_vstar")) { // Loop over relevant variable
		if (vars == "k1_vs_n") {
			k1 = step;
		} else if (vars == "ustar_vs_vstar") {
			vstar = step;
		}
		
		while ((k1 < k0 && vars == "k1_vs_n") || (vstar <= 5 && vars == "ustar_vs_vstar")) {

			if (vars == "k1_vs_n") {
				std::cout << " n: " << n << " K1: " << k1 << std::endl;
			} else if (vars == "ustar_vs_vstar") {
				std::cout << " ustar: " << ustar << " vstar: " << vstar << std::endl;
				k0 = 0.98 * 0.5 * (1/(ustar+vstar));
				k1 = 0.98 * k0;
			}	
			//Setting up the variables used in the DESolver
			desolver.set_k1(k1);
			desolver.set_n(n);
			desolver.set_ustar(ustar);
			desolver.set_vstar(vstar);
			desolver.set_k0(k0);
			desolver.set_alpha(1);
			desolver.set_y(1, 0);
			desolver.initialize();

			//Computing the fundemental matrix and the floquet multipliers
			fundemental_matrix = compute_fm(desolver);
			gsl_eigen_nonsymm(fundemental_matrix, floquet_multipliers, w);

			if (gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) > 1 || gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) > 1) {

				if (vars == "k1_vs_n") {
					out_file << n << " " << k1 << std::endl;
				} else if (vars == "ustar_vs_vstar") {
					out_file << ustar << " " << vstar << std::endl;
				}
				
				//Line below is for debugging
			/* std::cout << gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) << " " << gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) << std::endl; */
			}
			if (vars == "k1_vs_n") {
				k1 += step;
			} else if (vars == "ustar_vs_vstar") {
				vstar += step;
			}
		}
		if (vars == "k1_vs_n") {
			n += step;
		} else if (vars == "ustar_vs_vstar") {
			ustar += step;
		}
	}
}

void run_bifurcation(int argc, char* argv[], std::string vars) { // Still need to implement different vars

	//Check arguments are correct (might move this to its own function)
	if (vars == "alpha") {
		if (argc < 10) {
			printUsage(true);
			return;
		} else {
			//Parse arguments
			t0 = std::stoi(argv[2]);
			points = std::stoi(argv[3]);
			ustar = std::stod(argv[4]);
			vstar = std::stod(argv[5]);
			k0 = std::stod(argv[6]);
			k1 = std::stod(argv[7]);
			n = std::stod(argv[8]);
			step = std::stod(argv[9]);
		}
	} // Still need to implement this for different vars

	/*  else if (vars == "ustar_vs_vstar") { */
	/* 	if (argc < 5) { */
	/* 		printUsage(true); */
	/* 		return; */
	/* 	} else { */
	/* 		//Parse arguments */
	/* 		points = std::stoi(argv[2]); */
	/* 		n = std::stod(argv[3]); */
	/* 		step = std::stod(argv[4]); */
	/* 		ustar = step; */
	/* 	} */
	/* } */

	//Defining the linear DESolver (the parameters will be specificed later)
	alpha = 1;
	DESolver desolver(0);

	std::string command = "mkdir -p ../output/bifurcation/" + vars;
	std::system(command.c_str());
	std::ofstream out_file;
	if (vars == "alpha") {
		out_file.open("../output/bifurcation/" + vars + "/ustar=" + std::string(argv[4]) + "_vstar=" + std::string(argv[5]) + "_k0=" + std::string(argv[6]) + "_k1=" + std::string(argv[7]) + "_n=" + std::string(argv[8]) + ".dat");

	} //Still need to implement this for different vars
	/* } else if (vars == "ustar_vs_vstar") { */
/* 	out_file.open("../output/fm/" + vars + "/" + std::string(argv[3]) + ".dat"); */
/* } */

	while ( (alpha >= 0 && vars == "alpha") /*|| (ustar <= 5 && vars == "ustar_vs_vstar")*/) { // Loop over relevant variable

		if (vars == "alpha") {
			std::cout << "alpha=" << alpha << std::endl;

		} // Still need to implement this for different vars
		/* } else if (vars == "ustar_vs_vstar") { */
			/* 	std::cout << " ustar: " << ustar << " vstar: " << vstar << std::endl; */
			/* 	k0 = 0.98 * 0.5 * (1/(ustar+vstar)); */
			/* 	k1 = 0.98 * k0; */
		/* }	 */

		//Setting up the variables used in the DESolver
		desolver.set_k1(k1);
		desolver.set_n(n);
		desolver.set_ustar(ustar);
		desolver.set_vstar(vstar);
		desolver.set_k0(k0);
		desolver.set_alpha(alpha);
		desolver.set_y(ustar*(1+pow(10,-3)), vstar*(1+pow(10,-3)));
		desolver.initialize();
		tf = desolver.get_period();
		dt = tf/points;

		//Run the initial behavior
		desolver.solve(0, t0, dt);
		double current_u = desolver.get_y()[0];
		double current_v = desolver.get_y()[1];
		std::vector<std::pair<double, double>> points; // Points after transient
		points.clear(); // Just in case it was not cleared from the previous alpha value

		points.push_back(std::make_pair(current_u, current_v)); //This will be used to check if the initial time was not enough to capture the transient
		int repeat_index = 0;

		for (int i = 0; i < 100; i++) {

			desolver.solve(t0+i*tf, t0+(i+1)*tf, dt);
			current_u = desolver.get_y()[0];
			current_v = desolver.get_y()[1];

			// Check if the pattern started repeating after initial behavior
			repeat_index = 0;
			for (int j = 0; j < points.size(); j++) {
				if (std::abs(points[j].first - current_u) < 1e-6 && std::abs(points[j].second - current_v) < 1e-6) {
					repeat_index = j;
					break; // Repeating pattern found
				}
			}


			// If repetition is found, break the outer loop as well
			if (repeat_index != 0) break;
			
			// Add the current point to the list of points
			points.push_back({current_u, current_v});

		}

		for (int i = repeat_index; i < points.size(); i++) {
			out_file << alpha << " " << points[i].first << " " << points[i].second << std::endl;
		}

		if (vars == "alpha") {
			alpha -= step;
		} // Still need to implement this for different vars
		/* } else if (vars == "ustar_vs_vstar") { */
		/* 	vstar += step; */
		/* } */
	}
}


int main(int argc, char* argv[]) {
    if (argc > 1 && argv[1] != nullptr) {
		output_type = argv[1];
		if (output_type == "time_series") {
			run_time_series(argc, argv);
		} else if (output_type == "stability_k1_vs_n") {
			run_stability(argc,argv, "k1_vs_n");
		} else if (output_type == "stability_ustar_vs_vstar") {
			run_stability(argc,argv, "ustar_vs_vstar");
		} else if (output_type == "bifurcation_alpha") {
			run_bifurcation(argc, argv, "alpha");
			return 0;
		}
    } else {
        printUsage(true);
    }
}

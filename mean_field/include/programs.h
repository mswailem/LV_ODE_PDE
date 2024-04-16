#ifndef PROGRAMS_H
#define PROGRAMS_H

#include "DESolver.h"
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include "theory.h"
#include "iohandler.h"

// A struct to handle the logic of a varying parameter
struct VaryingParam {
	std::string name;
	double start;
	double end;
	double step;

	// Constructor
	VaryingParam(std::string name, double start, double end, double step) : name(name), start(start), end(end), step(step) {}

	// Get the number of points that this parameter will be varied over
	int get_num_of_points() {
		return ceil((end - start) / step);
	}
};

// Function to calculate the bifurcation diagram of the non-linear system for a given parameter
inline void bifurcation_diagram(std::unordered_map<std::string, double> p, VaryingParam v, double t0, int period_points) {
	
	DESolver desolver(0);

	int points = v.get_num_of_points();
	
	std::ofstream out_file = create_outfile("bifurcation", "bifurcation_test.dat");
	
	int progress = 0;
	std::cout << "\r"<< "Progress: " << progress << "/" << points << std::flush;

	for (int i = 0; i < points; i++) {
	
		desolver.set_params(p);
		desolver.set_y(p["us"]*(1+pow(10,-3)), p["vs"]*(1+pow(10,-3)));
		std::vector<std::pair<double, double>> fps = compute_stationary_points(desolver, t0, period_points, 1e-4);

		for (int j = 0; j < fps.size(); j++) {
			out_file << p[v.name] << " " << fps[j].first << " " << fps[j].second << std::endl;
		}

		p[v.name] += v.step;
		progress++;
		std::cout << "\r"<< "Progress: " << progress << "/" << points << std::flush;
	}
}

// Solve the non-linear system of ODEs
inline void time_series(std::unordered_map<std::string, double> p, double t0, double tf, double dt, std::vector<double> y0) {

	DESolver desolver(0);
	std::vector<double> y = y0;
	desolver.set_params(p);
	desolver.set_y(y[0], y[1]);
	desolver.initialize();

	std::ofstream out_file = create_outfile("time_series", "test.dat");

	for (double t = t0; t < tf; t += dt) {
		std::cout << "\r" << "Progress: " << t << "/" << tf << std::flush;
		out_file << t << " " << y[0] << " " << y[1] << std::endl;
		desolver.solve(t, t+dt, dt);
		y = desolver.get_y();
	}
	out_file << tf << " " << y[0] << " " << y[1] << std::endl;

}

// NOTE: Might have to change the implementation of this a little bit if the range of the v2 depends on the value of v1
// Function to calculate the stability diagram as two variables are varied
inline void stability(std::unordered_map<std::string, double> p, VaryingParam v1, VaryingParam v2, int period_points) {

	DESolver desolver(1);
	gsl_vector_complex *floquet_multipliers = gsl_vector_complex_alloc(2);
	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
	gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);
	std::ofstream out_file = create_outfile("stability", "stability_test.dat");

	int points1 = v1.get_num_of_points();
	int points2 = v2.get_num_of_points();
	int max_progress = points1 * points2;
	int progress = 0;

	std::cout << "\r" << "Progress: " << progress << "/" << max_progress << std::flush;
	for (int i = 0; i < points1; i++) {
		p[v2.name] = v2.start;
		for (int j = 0; j < points2; j++) {
			desolver.set_params(p);
			fundemental_matrix = compute_fm(desolver, period_points);
			gsl_eigen_nonsymm(fundemental_matrix, floquet_multipliers, w);
			bool unstable = gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) > 1 || gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) > 1;
			if (unstable) {
				out_file << p[v1.name] << " " << p[v2.name] << std::endl;
			}
			p[v2.name] += v2.step;
			progress++;
			std::cout << "\r" << "Progress: " << progress << "/" << max_progress << std::flush;
		}
		p[v1.name] += v1.step;
	}
	gsl_vector_complex_free(floquet_multipliers);
	gsl_matrix_free(fundemental_matrix);
	gsl_eigen_nonsymm_free(w);
}

inline void dispersion_relation(std::unordered_map<std::string, double> p, VaryingParam v) {

	std::ofstream out_file = create_outfile("dispersion", "test.dat");
	int points = v.get_num_of_points();
	int progress = 0;
	std::cout << "\r" << "Progress: " << progress << "/" << points << std::flush;

	for (int i = 0; i < points; i++) {
		gsl_vector_complex *eigenvalues = get_eigenvalues(p);

		gsl_complex lambda1 = gsl_vector_complex_get(eigenvalues, 0);
		gsl_complex lambda2 = gsl_vector_complex_get(eigenvalues, 1);
		out_file << p[v.name] << " " << GSL_REAL(lambda1) << " " << GSL_IMAG(lambda1) << " " << GSL_REAL(lambda2) << " " << GSL_IMAG(lambda2) << std::endl;
	}
}

#endif //PROGRAMS_H

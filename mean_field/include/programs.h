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
#include <omp.h>


// Function to calculate the bifurcation diagram of the non-linear system for a given parameter
inline void bifurcation_diagram(std::unordered_map<std::string, double> p, VaryingParam v, double t0, int period_points, std::string filename) {

	int points = v.get_num_of_points();
	
	std::ofstream out_file = create_outfile("bifurcation", filename);
	
	int progress = 0;
	std::cout << "\r"<< "Progress: " << progress << "/" << points << std::flush;

	#pragma omp parallel for shared(out_file, progress)
	for (int i = 0; i < points; i++) {

		std::unordered_map<std::string, double> p_local = p;
		double p_value = v.start + i*v.step;
		p_local[v.name] = p_value;

		DESolver desolver(0);

		desolver.set_params(p_local);
		desolver.set_y(p_local["us"]*(1+pow(10,-3)), p_local["vs"]*(1+pow(10,-3)));
		std::vector<std::pair<double, double>> fps = compute_stationary_points(desolver, t0, period_points, 1e-4);

		#pragma omp critical // Protects file writing 
		{
			for (int j = 0; j < fps.size(); j++) {
				out_file << p_value << " " << fps[j].first << " " << fps[j].second << std::endl;
			}
		}

		#pragma omp atomic
		++progress;

		std::cout << "\rProgress: " << progress << "/" << points << std::flush;
	}
}

// Creates a bifurcation phase plot for two varying parameters
inline void bifurcation_diagram(std::unordered_map<std::string, double> p, VaryingParam v1, VaryingParam v2, double t0, int period_points, std::string filename) {


	int points1 = v1.get_num_of_points();
	int points2 = v2.get_num_of_points();

	std::ofstream out_file = create_outfile("bifurcation", filename);

	int progress = 0;
	int max_progress = points1 * points2;
	std::cout << "\r"<< "Progress: " << progress << "/" << max_progress << std::flush;

	std::vector<std::pair<double, double>> fps;
	fps.clear(); // Just in case
	
	for (int i = 0; i < points1; i++) {
		#pragma omp parallel for shared(out_file, progress)
		for (int j = 0; j < points2; j++) {
			double p_value = v2.start + j*v2.step;
			std::unordered_map<std::string, double> p_local = p;
			p_local[v2.name] = p_value;
			DESolver desolver(0);
			desolver.set_params(p_local);
			desolver.set_y(p_local["us"]*(1+pow(10,-3)), p_local["vs"]*(1+pow(10,-3)));
			fps = compute_stationary_points(desolver, t0, period_points, 1e-4);

			int num_of_points = fps.size();
			#pragma omp critical // Protects file writing
			{
			out_file << p_local[v1.name] << " " << p_local[v2.name] << " " << num_of_points << std::endl;
			}

			#pragma omp atomic
			++progress;
			std::cout << "\r"<< "Progress: " << progress << "/" << max_progress << std::flush;
		}
		p[v1.name] += v1.step;
	}
}

// Solve the non-linear system of ODEs
inline void time_series(std::unordered_map<std::string, double> p, double t0, double tf, double dt, std::vector<double> y0, std::string filename) {

	DESolver desolver(0);
	std::vector<double> y = y0;
	desolver.set_params(p);
	desolver.set_y(y[0], y[1]);
	desolver.initialize();

	std::ofstream out_file = create_outfile("time_series", filename);

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
inline void stability(std::unordered_map<std::string, double> p, VaryingParam v1, VaryingParam v2, int period_points, std::string filename) {

	std::ofstream out_file = create_outfile("stability", filename);

	int points1 = v1.get_num_of_points();
	int points2 = v2.get_num_of_points();
	int max_progress = points1 * points2;
	int progress = 0;

	std::cout << "\r" << "Progress: " << progress << "/" << max_progress << std::flush;
	for (int i = 0; i < points1; i++) {

		#pragma omp parallel for shared(out_file, progress)
		for (int j = 0; j < points2; j++) {

			gsl_vector_complex *floquet_multipliers = gsl_vector_complex_alloc(2);
			gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
			gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);

			DESolver desolver(1);
			double p_value = v2.start + j*v2.step;
			std::unordered_map<std::string, double> p_local = p;
			p_local[v2.name] = p_value;

			desolver.set_params(p_local);
			compute_fm(desolver, period_points, fundemental_matrix);
			gsl_eigen_nonsymm(fundemental_matrix, floquet_multipliers, w);
			bool unstable = gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) > 1 || gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) > 1;

			if (unstable) {
				#pragma omp critical // Protects file writing 
				{
				out_file << p_local[v1.name] << " " << p_local[v2.name] << std::endl;
				}
			}

			#pragma omp atomic
			++progress;
			std::cout << "\r" << "Progress: " << progress << "/" << max_progress << std::flush;
			gsl_vector_complex_free(floquet_multipliers);
			gsl_matrix_free(fundemental_matrix);
			gsl_eigen_nonsymm_free(w);
		}
		p[v1.name] += v1.step;
	}
}

inline void dispersion_relation(std::unordered_map<std::string, double> p, VaryingParam v, std::string filename) {

	std::ofstream out_file = create_outfile("dispersion", filename);
	int points = v.get_num_of_points();
	int progress = 0;
	std::cout << "\r" << "Progress: " << progress << "/" << points << std::flush;

	#pragma omp parallel for shared(out_file, progress)
	for (int i = 0; i < points; i++) {

		double p_value = v.start + i*v.step;
		std::unordered_map<std::string, double> p_local = p;
		p_local[v.name] = p_value;

		gsl_vector_complex *eigenvalues = get_eigenvalues(p_local);

		gsl_complex lambda1 = gsl_vector_complex_get(eigenvalues, 0);
		gsl_complex lambda2 = gsl_vector_complex_get(eigenvalues, 1);

		#pragma omp critical // Protects file writing 
		{
			out_file << p_local[v.name] << " " << GSL_REAL(lambda1) << " " << GSL_IMAG(lambda1) << " " << GSL_REAL(lambda2) << " " << GSL_IMAG(lambda2) << std::endl;
		}
		#pragma omp atomic
		++progress;
		std::cout << "\r" << "Progress: " << progress << "/" << points << std::flush;
		gsl_vector_complex_free(eigenvalues);
	}
}

#endif //PROGRAMS_H

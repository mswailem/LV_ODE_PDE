#ifndef THEORY_H
#define THEORY_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex_double.h>
#include <unordered_map>
#include "DESolver.h"

// Function to compute the long-time stationary points of the system
inline std::vector<std::pair<double, double>> compute_stationary_points(DESolver &solver, double t0, int points, double tolerance) {
	solver.initialize();
	double tf = solver.get_period();
	double dt = tf / points;
	solver.solve(0, t0, dt);
	double current_u = solver.get_y()[0];
	double current_v = solver.get_y()[1];
	std::vector<std::pair<double, double>> fps; // Points after transient
	fps.clear();

	fps.push_back(std::make_pair(current_u, current_v)); //This will be used to check if the initial time was not enough to capture the transient
	int repeat_index = 0;

	for (int i = 0; i < 8; i++) {

		solver.solve(t0 + i * tf, t0+(i+1)*tf, dt); // Solve for one period
		
		current_u = solver.get_y()[0];
		current_v = solver.get_y()[1];

		// Check if the pattern started repeating after initial behavior
		repeat_index = -1;
		for (int j = 0; j < fps.size(); j++) {
			if (std::abs(fps[j].first - current_u) < tolerance && std::abs(fps[j].second - current_v) < tolerance) {
				repeat_index = j;
				break; // Repeating pattern found
			}
		}

		// If repetition is found, break the outer loop as well
		if (repeat_index != -1) {
			// Erase all points before the repeating pattern, keeping the repeating pattern intact
			fps.erase(fps.begin(), fps.begin() + repeat_index);
			break; // Exit the loop as we've found the repeating pattern
		}
		// Add the current point to the list of points
		fps.push_back({current_u, current_v});

	}
	return fps;
}

// Compute the eigenvalues of the Jacobian matrix
inline gsl_vector_complex* get_eigenvalues(std::unordered_map<std::string, double> p) {

	gsl_matrix *Jacobian = gsl_matrix_alloc(2, 2);
	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
	
	double a = -p["du"]*p["wn"]*p["wn"];
	double b = -(p["us"]+p["vs"])*p["k0"]+1;
	double c = (p["vs"]*p["vs"]/p["us"])*p["k0"]-(p["vs"]/p["us"]);
	double d = -p["dv"]*p["wn"]*p["wn"]-p["vs"]*p["k0"];

	gsl_matrix_set(Jacobian, 0, 0, a);
	gsl_matrix_set(Jacobian, 0, 1, b);
	gsl_matrix_set(Jacobian, 1, 0, c);
	gsl_matrix_set(Jacobian, 1, 1, d);

	gsl_vector_complex *eigenvalues = gsl_vector_complex_alloc(2);
	
	gsl_eigen_nonsymm(Jacobian, eigenvalues, w);

	gsl_matrix_free(Jacobian);
	gsl_eigen_nonsymm_free(w);
	return eigenvalues;
}

inline gsl_matrix* compute_fm(DESolver& solver, double points) { 
		// Variable defintions
		double tf = solver.get_period();
		double dt = tf/points;

		//Allocating memory for gsl variables
		static gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);

		//Solving the system for the the (1,0) vector
		solver.set_y(1, 0);
		solver.initialize();
		solver.solve(0, tf, dt);
		std::vector<double> y = solver.get_y();

		//Assigning the values of the fundemental matrix
		gsl_matrix_set(fundemental_matrix, 0, 0, y[0]);
		gsl_matrix_set(fundemental_matrix, 0, 1, y[1]);

		//Solving the system for the (0,1) vector
		solver.set_y(0, 1);
		solver.initialize();
		solver.solve(0, tf, dt);
		y = solver.get_y();

		//Assigning the values of the fundemental matrix
		gsl_matrix_set(fundemental_matrix, 1, 0, y[0]);
		gsl_matrix_set(fundemental_matrix, 1, 1, y[1]);

		return fundemental_matrix;
}

#endif //THEORY_H

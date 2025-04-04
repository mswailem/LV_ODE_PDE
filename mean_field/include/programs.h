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


// TODO: I added the option to handle the cases where I want k0 and k1 to run to their maximum values, check what I have done in the stability program, also implement this in the other programs


// Function to calculate the bifurcation diagram of the non-linear system for a given parameter
inline void bifurcation_diagram(std::unordered_map<std::string, double> p, VaryingParam v, double t0, std::string filename) {

	std::ofstream out_file = create_outfile("bifurcation", filename);
	
	VaryingParam v_local = handle_max_inputs(p, v);

	std::vector<double> y0 = {p["us"]+p["da"], 1+p["db"]};

	for (auto &v_value : v_local.values) {
		std::unordered_map<std::string, double> p_local = p;
		p_local[v.name] = v_value;

		DESolver desolver(0);

		desolver.set_params(p_local);
		if (!desolver.check_param_region()) {
			continue;
		}
		desolver.set_y(y0[0], y0[1]);
		std::vector<std::pair<double, double>> fps = compute_stationary_points(desolver, 1e-3);
		y0 = {fps[0].first + p_local["da"], fps[0].second + p_local["db"]};

		for (int j = 0; j < fps.size(); j++) {
			out_file << v_value << " " << fps[j].first << " " << fps[j].second << std::endl;
		}
	}
}

// Creates a bifurcation phase plot for two varying parameters
inline void bifurcation_diagram(std::unordered_map<std::string, double> p, VaryingParam v1, VaryingParam v2, double t0, std::string filename) {

    std::vector<std::string> results; // Vector to store results for file output
	VaryingParam v1_local = handle_max_inputs(p, v1);
	std::vector<double> y0;


    for (auto &v1_value : v1_local.values) {
		p[v1.name] = v1_value;
        std::vector<std::string> local_results; // Local vector for each outer loop iteration
		VaryingParam v2_local = handle_max_inputs(p, v2);
		y0 = {p["us"]+p["da"], 1+p["db"]};

        for (auto &v2_value : v2_local.values) {
			std::cout << "\r" << "Progress: " << v1_value << " " << v2_value << "/" << v1_local.end << " " << v2_local.end << std::flush;
            std::unordered_map<std::string, double> p_local = p;
            p_local[v1.name] = v1_value;
            p_local[v2.name] = v2_value;

            DESolver desolver(0);
            desolver.set_params(p_local);
			if (!desolver.check_param_region()) {
				continue;
			}
            desolver.set_y(y0[0], y0[1]);
            std::vector<std::pair<double, double>> fps = compute_stationary_points(desolver, 1e-3);
			y0 = {fps[0].first + p_local["da"], fps[0].second + p_local["db"]};

            std::string result = std::to_string(p_local[v1.name]) + " " + std::to_string(p_local[v2.name]) + " " + std::to_string(fps.size()) + "\n";
			
            local_results.push_back(result);
        }

        // Combine local results into the main results vector
        results.insert(results.end(), local_results.begin(), local_results.end());
    }

    // Write all results to file sequentially outside the parallel region
    std::ofstream out_file = create_outfile("bifurcation", filename);
    for (const auto& line : results) {
        out_file << line;
    }
    out_file.close();
}

// Solve the non-linear system of ODEs
inline void time_series(std::unordered_map<std::string, double> p, double t0, double tf, std::vector<double> y0, std::string filename) { // dt is the save frequency of my program, not the timestep used to solve the system

	DESolver desolver(0);
	std::vector<double> y = y0;
	if (p["k1"] < 0) { 
		p["k1"] = 0.98*max_k1(p["us"], p["k0"], p["alpha"]);
	}
	desolver.set_params(p);
	if (!desolver.check_param_region()) {
		throw std::runtime_error("Parameters are outside of the physical oscillatoy region");
	}
	desolver.set_y(y[0], y[1]);
	desolver.initialize();
	double dt = desolver.get_period()/100;
	int num_of_periods = ceil(tf);

	std::ofstream out_file = create_outfile("time_series", filename);

	for (int i = 0; i < num_of_periods*100; i++) {
		double t = t0 + i*dt;
		std::cout << "\r" << "Progress: " << i << "/" << num_of_periods*100 << std::flush;
		out_file << t << " " << y[0] << " " << y[1] << std::endl;
		desolver.solve(t, t+dt);
		y = desolver.get_y();
	}
	out_file << tf << " " << y[0] << " " << y[1] << std::endl;

}

// NOTE: Might have to change the implementation of this a little bit if the range of the v2 depends on the value of v1
// Function to calculate the stability diagram as two variables are varied
inline void stability(std::unordered_map<std::string, double> p, VaryingParam v1, VaryingParam v2, std::string filename) {
	std::ofstream out_file = create_outfile("stability", filename);

	// Stability analysis only makes sense for alpha = 1
	p["alpha"] = 1;

	VaryingParam v1_local = handle_max_inputs(p, v1);

	for (auto &v1_value : v1_local.values) {
		p[v1.name] = v1_value;
		VaryingParam v2_local = handle_max_inputs(p, v2);

		/* #pragma omp parallel for shared(out_file, progress) */

		#pragma omp parallel for shared(out_file)
		for (auto &v2_value : v2_local.values) {
			std::cout << "\r" << "Progress: " << v1_value << "/" << v1_local.end << ", " << v2_value << "/" << v2_local.end << std::flush;
			gsl_vector_complex *floquet_multipliers = gsl_vector_complex_alloc(2);
			gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
			gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);

			DESolver desolver(1);
			std::unordered_map<std::string, double> p_local = p;
			p_local[v2.name] = v2_value;

			// TODO: Find a better way to do this, maybe I can define a function that does somewhere and then call it here?
			// Check if k1 should be equal to close maximum value or not
			
			if (p_local["k1"] < 0) { 
				p_local["k1"] = 0.98*max_k1(p_local["us"], p_local["k0"], p_local["alpha"]);
			}

			desolver.set_params(p_local);
			if (!desolver.check_param_region()) {
				gsl_vector_complex_free(floquet_multipliers);
				gsl_matrix_free(fundemental_matrix);
				gsl_eigen_nonsymm_free(w);
				continue;
			}
			compute_fm(desolver, fundemental_matrix);
			gsl_eigen_nonsymm(fundemental_matrix, floquet_multipliers, w);
			bool unstable = gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) > 1 || gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) > 1;

			if (unstable) {
				#pragma omp critical // Protects file writing 
				{
				out_file << p_local[v1.name] << " " << p_local[v2.name] << std::endl;
				}
			}

			gsl_vector_complex_free(floquet_multipliers);
			gsl_matrix_free(fundemental_matrix);
			gsl_eigen_nonsymm_free(w);
		}
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

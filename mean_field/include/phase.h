#ifndef PHASE_H
#define PHASE_H

#include "DESolver.h"
#include <vector>

inline std::vector<std::pair<double, double>> compute_fixed_points(DESolver &solver, double t0, int points, double tolerance) {
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

	for (int i = 0; i < 64; i++) {

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

#endif //PHASE_H

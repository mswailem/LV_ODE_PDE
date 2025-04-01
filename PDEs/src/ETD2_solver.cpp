#include "ETD2_solver.h"
#include <cmath>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Constructor function, set member variables and initialize FFTW arrays and plans
ETD2_solver::ETD2_solver(double t_start_init, int t_end_init, double L_init, int N_init, double ustar_init, double k0_init, double k1_init, double n_init, double alpha_init , double D_init, double tt_spacing)

	// This is a "member initializer list" to initialize the member variables
    : t_start(t_start_init), num_of_periods(t_end_init), L(L_init), N(N_init), astar(ustar_init), k0(k0_init), k1(k1_init), n(n_init), alpha(alpha_init), D(D_init), t_spacing(tt_spacing),
	  c_a(N * N), c_b(N * N), f1(N * N),
      g1(N * N), f2(N * N), g2(N * N), f3(N * N), g3(N * N), f4(N * N),
      g4(N * N), xx(N), yy(N), kx(N), ky(N) {

	// Calculate forcing period
	omega0 = calculate_omega0();
	double forcing_omega = omega0/n;
	forcing_period = 2 * M_PI / forcing_omega;

	// Calculate time and space spacing related variables
	t_spacing = t_spacing * forcing_period; // User specifies t_spacing in terms of a multiple of the forcing period (I think this will work even if it is not an integer multiple)
	dt = 0.005 * forcing_period; // Time stepping is set to 0.005 of the forcing period to ensure enough smapling of one forcing period
	dx = L / N;
	t_end = (double)num_of_periods * forcing_period;
	t_points = std::ceil((t_end - t_start) / dt);

	// Initialize arrays for FFTW
	a = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	b = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	f = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	g = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	f_prev = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	g_prev = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);

	// Initialize FFT plans
	forward_transform = fftw_plan_dft_2d(N, N, a, a, FFTW_FORWARD, FFTW_MEASURE);
	backward_transform =
		fftw_plan_dft_2d(N, N, a, a, FFTW_BACKWARD, FFTW_MEASURE);
}

// Deconstructor function, clean up FFTW arrays and plans
ETD2_solver::~ETD2_solver() {
  fftw_destroy_plan(forward_transform);
  fftw_destroy_plan(backward_transform);
  fftw_free(a);
  fftw_free(b);
  fftw_free(f);
  fftw_free(g);
  fftw_free(f_prev);
  fftw_free(g_prev);
}

// Calculate omega0
double ETD2_solver::calculate_omega0() {
  double term1 = (1 / astar) * (k0 - 1) * (k0 - 1);
  double term2 = 0.25 * k0 * (-4 + 3 * k0);
  return sqrt(term1 + term2);
}

// Calculate the predation rate
double ETD2_solver::lambda(double t) {
	double lambda1 = (1 / astar) * (1 - k(t)) - k(t);
	double lambda0 = (1 / astar) * (1 - k0) - k0;
	return (1 - alpha) * lambda0 + alpha * lambda1;
}

// Calculate the forcing function
double ETD2_solver::k(double t) {
	if (n == 0) { // I specify n = 0 to indicate a constant environment
		return k0;
	}
	return k0 + k1 * cos((1 / n) * t * omega0);
}

// Calculate the wavenumbers for FFT
void ETD2_solver::calculate_wavenumbers() {
  double dk = 2 * M_PI / L; // Wavenumber spacing
  int N_half = N / 2;

  for (int i = 0; i < N; i++) {
    // Adjust for negative frequencies
    kx[i] = (i < N_half) ? i * dk : (i - N) * dk;
    ky[i] = (i < N_half) ? i * dk : (i - N) * dk;
  }
}

// Calculates the ETD2 coefficients for the ETD2 method in log space
void ETD2_solver::calculate_coeffcients_log() {
	// These are temporary variables to hold -D*k^2 for each grid point in k-space
	double k2_a;
	double k2_b;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			k2_a = -D * kx[i] * kx[i] - D * ky[j] * ky[j];
			k2_b = -D * kx[i] * kx[i] - D * ky[j] * ky[j];

			// These are the standard formulas of the ETD2 method, with checks to avoid division by zero when c_a or c_b are zero, then the coefficients are set to their limit as c_a/c_b -> 0
			c_a[i * N + j] = k2_a;
			f1[i * N + j] = exp(c_a[i * N + j] * dt);
			f2[i * N + j] = c_a[i * N + j] == 0 ? dt : (f1[i * N + j] - 1) / c_a[i * N + j];
			f3[i * N + j] = c_a[i * N + j] == 0 ? 1.5 * dt : ((1 + c_a[i * N + j] * dt) * f1[i * N + j] - 1 - 2 * c_a[i * N + j] * dt) / (std::pow(c_a[i * N + j], 2) * dt);
			f4[i * N + j] = c_a[i * N + j] == 0 ? -0.5 * dt : (-f1[i * N + j] + 1 + c_a[i * N + j] * dt) / (std::pow(c_a[i * N + j], 2) * dt);

			c_b[i * N + j] = k2_b;
			g1[i * N + j] = exp(c_b[i * N + j] * dt);
			g2[i * N + j] = c_b[i * N + j] == 0 ? dt : (g1[i * N + j] - 1) / c_b[i * N + j];
			g3[i * N + j] = c_b[i * N + j] == 0 ? 1.5 * dt : ((1 + c_b[i * N + j] * dt) * g1[i * N + j] - 1 - 2 * c_b[i * N + j] * dt) / (std::pow(c_b[i * N + j], 2) * dt);
			g4[i * N + j] = c_b[i * N + j] == 0 ? -0.5 * dt : (-g1[i * N + j] + 1 + c_b[i * N + j] * dt) / (std::pow(c_b[i * N + j], 2) * dt);
		}
	}
}

// Calculates the ETD2 coefficients for the ETD2 method in log space (outdated method)
void ETD2_solver::calculate_coeffcients_regular() { 
	// These are temporary variables to hold -D*k^2 for each grid point in k-space
	double k2_u;
	double k2_v;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			k2_u = -D * kx[i] * kx[i] - D * ky[j] * ky[j];
			k2_v = -D * kx[i] * kx[i] - D * ky[j] * ky[j];

			// These are the standard formulas of the ETD2 method, with checks to avoid division by zero when c_a or c_b are zero, then the coefficients are set to their limit as c_a/c_b -> 0
			c_a[i * N + j] = k2_u - lambda(0);
			f1[i * N + j] = exp(c_a[i * N + j] * dt);
			f2[i * N + j] = c_a[i * N + j] == 0 ? dt : (f1[i * N + j] - 1) / c_a[i * N + j];
			f3[i * N + j] = c_a[i * N + j] == 0 ? 1.5 * dt : ((1 + c_a[i * N + j] * dt) * f1[i * N + j] - 1 - 2 * c_a[i * N + j] * dt) / (std::pow(c_a[i * N + j], 2) * dt);
			f4[i * N + j] = c_a[i * N + j] == 0 ? -0.5 * dt : (-f1[i * N + j] + 1 + c_a[i * N + j] * dt) / (std::pow(c_a[i * N + j], 2) * dt);

			c_b[i * N + j] = k2_v + 1;
			g1[i * N + j] = exp(c_b[i * N + j] * dt);
			g2[i * N + j] = c_b[i * N + j] == 0 ? dt : (g1[i * N + j] - 1) / c_b[i * N + j];
			g3[i * N + j] = c_b[i * N + j] == 0 ? 1.5 * dt : ((1 + c_b[i * N + j] * dt) * g1[i * N + j] - 1 - 2 * c_b[i * N + j] * dt) / (std::pow(c_b[i * N + j], 2) * dt);
			g4[i * N + j] = c_b[i * N + j] == 0 ? -0.5 * dt : (-g1[i * N + j] + 1 + c_b[i * N + j] * dt) / (std::pow(c_b[i * N + j], 2) * dt);
		}
	}
}

// Computes the magnitude of a complex number
double ETD2_solver::complex_magnitude(fftw_complex z) {
  return sqrt(z[0] * z[0] + z[1] * z[1]);
}

// Multiplies two complex numbers and saves the result in the third argument
void ETD2_solver::complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result) {
  fftw_complex result_temp;
  result_temp[0] = z[0] * w[0] - z[1] * w[1];
  result_temp[1] = z[0] * w[1] + z[1] * w[0];
  result[0] = result_temp[0];
  result[1] = result_temp[1];
}

// Calculates the non-linear terms in log space
void ETD2_solver::calculate_non_linear_log(double t) {
	// The non-linear terms include spatial derivatives, and so these variables will hold the first order spatial deivatives
	fftw_complex *deriv_a_x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex *deriv_b_x = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex *deriv_a_y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex *deriv_b_y = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);

	// Transform the system to k-space to compute spatial derivatives
	transform_to_k_space();

	// Compute the spatial derivatives in k-space
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			deriv_a_x[i * N + j][0] = kx[i] * a[i * N + j][1];
			deriv_a_x[i * N + j][1] = -kx[i] * a[i * N + j][0];
			deriv_b_x[i * N + j][0] = kx[i] * b[i * N + j][1];
			deriv_b_x[i * N + j][1] = -kx[i] * b[i * N + j][0];
			deriv_a_y[i * N + j][0] = ky[j] * a[i * N + j][1];
			deriv_a_y[i * N + j][1] = -ky[j] * a[i * N + j][0];
			deriv_b_y[i * N + j][0] = ky[j] * b[i * N + j][1];
			deriv_b_y[i * N + j][1] = -ky[j] * b[i * N + j][0];
		}
	}

	// Transform back to x-space
	transform_to_x_space();
	fftw_execute_dft(backward_transform, deriv_a_x, deriv_a_x);
	fftw_execute_dft(backward_transform, deriv_b_x, deriv_b_x);
	fftw_execute_dft(backward_transform, deriv_a_y, deriv_a_y);
	fftw_execute_dft(backward_transform, deriv_b_y, deriv_b_y);

	// Compute the non-linear terms
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f[i * N + j][0] = D * std::pow(deriv_a_x[i * N + j][0], 2) + D * std::pow(deriv_a_y[i * N + j][0], 2) + lambda(t) * exp(b[i * N + j][0]) - lambda(t);
			g[i * N + j][0] = D * std::pow(deriv_b_x[i * N + j][0], 2) + D * std::pow(deriv_b_y[i * N + j][0], 2) - lambda(t) * exp(a[i * N + j][0]) + (1 - k(t) * exp(a[i * N + j][0]) - k(t) * exp(b[i * N + j][0]));
			f[i * N + j][1] = 0;
			g[i * N + j][1] = 0;
		}
	}

	// Free derivative arrays
	fftw_free(deriv_a_x);
	fftw_free(deriv_b_x);
	fftw_free(deriv_a_y);
	fftw_free(deriv_b_y);
}

// Calculates the non-linear terms in regular space (outdated method)
void ETD2_solver::calculate_non_linear_regular() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      f[i * N + j][0] = lambda(0) * b[i * N + j][0] * a[i * N + j][0];
      g[i * N + j][0] = -lambda(0) * b[i * N + j][0] * a[i * N + j][0] - (a[i * N + j][0] * b[i * N + j][0] + std::pow(b[i * N + j][0], 2));
      f[i * N + j][1] = 0;
      g[i * N + j][1] = 0;
    }
  }
}

// Sets the initial conditions for the system based on the specified type
void ETD2_solver::set_initial_conditions(std::string type) {
	// This is only used for random initial conditions
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> d(0, 1e-1); // The second argument here specifies the size of the initial spatial fluctuations

	for (int i = 0; i < N; i++) {
		xx[i] = i * dx;
		for (int j = 0; j < N; j++) {
			if (i == 0) {
				yy[j] = j * dx;
			}
			if (type == "checkboard") { // Checkboard pattern
				if (xx[i] < (double)L / 2 && yy[j] < (double)L / 2) {
					a[i * N + j][0] = 1;
					a[i * N + j][1] = 0;
					b[i * N + j][0] = 0.01;
					b[i * N + j][1] = 0;
				} else if (xx[i] < (double)L / 2 && yy[j] >= (double)L / 2) {
					a[i * N + j][0] = 0.01;
					a[i * N + j][1] = 0;
					b[i * N + j][0] = 1;
					b[i * N + j][1] = 0;
				} else if (xx[i] >= (double)L / 2 && yy[j] < (double)L / 2) {
					a[i * N + j][0] = 0.01;
					a[i * N + j][1] = 0;
					b[i * N + j][0] = 1;
					b[i * N + j][1] = 0;
				} else {
					a[i * N + j][0] = 1;
					a[i * N + j][1] = 0;
					b[i * N + j][0] = 0.01;
					b[i * N + j][1] = 0;
				}
			} else if (type == "continuous") { // Continuous initial conditions with sines and cosines
				a[i * N + j][0] = exp(sin(2 * M_PI * xx[i] / L));
				a[i * N + j][1] = 0;
				b[i * N + j][0] = exp(cos(2 * M_PI * yy[j] / L));
				b[i * N + j][1] = 0;
			} else if (type == "constant") { // Constant initial conditions
				a[i * N + j][0] = astar + 0.01;
				a[i * N + j][1] = 0;
				b[i * N + j][0] = 1 + 0.01;
				b[i * N + j][1] = 0;
			} else if (type == "random") { // Staring from the fixed point of the average system with random perturbations
				double rand1 = d(gen);
				double rand2 = d(gen);
				if (rand1 < 0) {
					rand1 = -rand1;
				}
				if (rand2 < 0) {
					rand2 = -rand2;
				}
				a[i * N + j][0] = astar + rand1;
				a[i * N + j][1] = 0;
				b[i * N + j][0] = 1 + rand2;
				b[i * N + j][1] = 0;
			}
		}
	}
}

// Writes the data to a file
void ETD2_solver::write_data(std::string file_name, bool k_space) {

	// Output on k-space
	if (k_space) {
		std::ofstream file("../output/" + file_name + "_k.dat");
		file << std::scientific << std::setprecision(14); // Add more digits, 14 digits is the maximum precision for double in C++
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				file << kx[i] << " " << ky[j] << " " << complex_magnitude(a[i * N + j]) << " " << complex_magnitude(b[i * N + j]) << std::endl; // Data is outputted as a line: kx ky a(kx,ky) b(kx,ky)
			}
			file << std::endl;
		}
		file.close();
	} else { // Else output in x-space
		std::ofstream file("../output/" + file_name + "_x.dat");
		file << std::scientific << std::setprecision(14); // Add more digits, 14 digits is the maximum precision for double in C++
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				file << xx[i] << " " << yy[j] << " " << a[i * N + j][0] << " " << b[i * N + j][0] << std::endl; // Data is outputted as a line: x y a(x,y) b(x,y)
			}
			file << std::endl;
		}
		file.close();
	}
}

// Transform to k-space using FFTW
void ETD2_solver::transform_to_k_space() {
	fftw_execute_dft(forward_transform, a, a);
	fftw_execute_dft(forward_transform, b, b);
	fftw_execute_dft(forward_transform, f, f);
	fftw_execute_dft(forward_transform, g, g);

	// Normalize after FFT
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i * N + j][0] /= N * N;
			a[i * N + j][1] /= N * N;
			b[i * N + j][0] /= N * N;
			b[i * N + j][1] /= N * N;
			f[i * N + j][0] /= N * N;
			g[i * N + j][0] /= N * N;
			f[i * N + j][1] /= N * N;
			g[i * N + j][1] /= N * N;
		}
	}
}

// Transform back to x-space using FFTW
void ETD2_solver::transform_to_x_space() {
  fftw_execute_dft(backward_transform, a, a);
  fftw_execute_dft(backward_transform, b, b);
  fftw_execute_dft(backward_transform, f, f);
  fftw_execute_dft(backward_transform, g, g);
}

// Transform the system to log space
void ETD2_solver::transform_to_log_space() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      a[i * N + j][0] = log(a[i * N + j][0]);
      a[i * N + j][1] = 0;
      b[i * N + j][0] = log(b[i * N + j][0]);
      b[i * N + j][1] = 0;
    }
  }
}

// Transform the system back to regular space
void ETD2_solver::transform_to_regular_space() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      a[i * N + j][0] = exp(a[i * N + j][0]);
      a[i * N + j][1] = 0;
      b[i * N + j][0] = exp(b[i * N + j][0]);
      b[i * N + j][1] = 0;
    }
  }
}

// This is the update function for the time-stepping in the ETD2 method
void ETD2_solver::time_step(bool first_time) {
	// If this is the first iteration use ETD1
	if (first_time) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i * N + j][0] = f1[i * N + j] * a[i * N + j][0] + f2[i * N + j] * f[i * N + j][0];
				a[i * N + j][1] = f1[i * N + j] * a[i * N + j][1] + f2[i * N + j] * f[i * N + j][1];
				b[i * N + j][0] = g1[i * N + j] * b[i * N + j][0] + g2[i * N + j] * g[i * N + j][0];
				b[i * N + j][1] = g1[i * N + j] * b[i * N + j][1] + g2[i * N + j] * g[i * N + j][1];
				f_prev[i * N + j][0] = f[i * N + j][0];
				f_prev[i * N + j][1] = f[i * N + j][1];
				g_prev[i * N + j][0] = g[i * N + j][0];
				g_prev[i * N + j][1] = g[i * N + j][1];
			}
		}
	} else {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i * N + j][0] = f1[i * N + j] * a[i * N + j][0] + f3[i * N + j] * f[i * N + j][0] + f4[i * N + j] * f_prev[i * N + j][0];
				a[i * N + j][1] = f1[i * N + j] * a[i * N + j][1] + f3[i * N + j] * f[i * N + j][1] + f4[i * N + j] * f_prev[i * N + j][1];
				b[i * N + j][0] = g1[i * N + j] * b[i * N + j][0] + g3[i * N + j] * g[i * N + j][0] + g4[i * N + j] * g_prev[i * N + j][0];
				b[i * N + j][1] = g1[i * N + j] * b[i * N + j][1] + g3[i * N + j] * g[i * N + j][1] + g4[i * N + j] * g_prev[i * N + j][1];
				f_prev[i * N + j][0] = f[i * N + j][0];
				f_prev[i * N + j][1] = f[i * N + j][1];
				g_prev[i * N + j][0] = g[i * N + j][0];
				g_prev[i * N + j][1] = g[i * N + j][1];
			}
		}
	}
}

// Solve the system using log space
void ETD2_solver::solve_in_log() {
	double t = t_start; // Initialize time

	// Calculate variables that can be precomputed
	calculate_wavenumbers();
	calculate_coeffcients_log();

	// Set initial conditions and write data at t=0
	set_initial_conditions(std::string("random"));
	write_data("0");
	transform_to_k_space();
	write_data("0", true);
	transform_to_x_space();

	// Calculate the non-linear terms in log space
	transform_to_log_space();
	calculate_non_linear_log(t);
	transform_to_k_space(); // Time update needs to happen in k-space
	
	int counter = 1; // This counter is used to know when it is time to log the data to a file
	for (int i = 0; i < t_points; i++) {
		t += dt;
		time_step(i == 0);
		transform_to_x_space(); // Transform back to x-space to calculate the non-linear terms/writing to file
		if (t >= (double)counter * t_spacing) {
			std::cout << "Writing data at t = " << t << "/" << t_end << "\r" << std::flush;
			transform_to_regular_space();
			write_data(std::to_string(counter));
			transform_to_k_space();
			write_data(std::to_string(counter), true);
			transform_to_x_space();
			transform_to_log_space();
			counter += 1;
		}
		calculate_non_linear_log(t);
		transform_to_k_space();
	}
}

// Solve the system in regular space (outdated method, not recommended to use)
void ETD2_solver::solve_in_regular() {
	double t = t_start;
	calculate_wavenumbers();
	calculate_coeffcients_regular();
	set_initial_conditions(std::string("constant"));
	write_data("0");
	calculate_non_linear_regular();
	transform_to_k_space();
	/* write_data("0", true); */
	int counter = 1;
	for (int i = 0; i < t_points; i++) {
		t += dt;
		time_step(i == 0);
		transform_to_x_space();
		if (t >= counter) {
			write_data(std::to_string(counter));
			counter++;
		}
		calculate_non_linear_regular();
		transform_to_k_space();
		/* if (t >= counter) { */
		/* 	write_data(std::to_string(counter), true); */
		/* 	counter++; */
		/* } */
	}
}

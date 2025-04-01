#ifndef ETD2_H
#define ETD2_H

#include <fftw3.h>
#include <string>
#include <vector>
#include <random>

// This class holds all the necessary methods and variables to solve partial differential equations using the Exponential time differencing to the second order (ETD2) method.
class ETD2_solver {
	public:
		ETD2_solver(double t_start_init, int t_end_init, double L_init, int N_init, double ustar_init, double k0_init, double k1_init, double n_init, double alpha_init, double D_init, double tt_spacing); // Constructor
		double complex_magnitude(fftw_complex z); // Computes the magnitude of a complex number
		void complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result); // Multiplies two complex numbers
		void solve_in_log(); // Solves the PDE in logarithmic space (the results are produced in regular space)
		void solve_in_regular(); // Solves the PDE in regular space (outdated method, and I haven't really checked it in a while)
		~ETD2_solver();

	private:
	// Precomputing terms
	void calculate_wavenumbers(); // Precomputes the wavenumbers for the Fourier space
	void calculate_coeffcients_log(); // Precomputes the coefficients of the ETD2 method when equations are implemented in log space
	void calculate_coeffcients_regular(); // Precomputes the coefficients of the ETD2 method when equations are implemented in regular space (this method is outdated and hasn't been checked in a while)
	void calculate_non_linear_log(double t); // Calculates the non-linear terms for the equation in log space
	void calculate_non_linear_regular(); // Calculates the non-linear terms for the equation in regular space (this method is outdated and hasn't been checked in a while)
	double calculate_omega0(); // Calculates the intrinsic frequence omega0

	// Time stepping methods
	void set_initial_conditions(std::string type); // Sets initial conditions, available types: "checkboard" (checkerboard pattern), "continuous" (sine/cosine), "constant" (constant values), "random" (random perturbations around the average environment fixed point)
	void time_step(bool first_time); // This is the update function that actually implements the ETD2 equations, the bool "first_time" is used since ETD2 first time step is implemented as ETD1 since ETD2 depends on two previous time steps.

	// Transformation methods
	void transform_to_k_space(); // Transform the system to Fourier space (k-space) using FFTW
	void transform_to_x_space(); // Transform the system back to real space (x-space) using FFTW
	void transform_to_log_space(); // Transform the system to log space
	void transform_to_regular_space(); // Transform the system back to regular space

	void write_data(std::string file_name, bool k_space = false); // Writes the data to a file, if k_space is true it writes the data in Fourier space.

	//Variables
	
	//Time variables
	double t_start; // Starting time
	double t_end; // Final time
	int num_of_periods; // The number of forcing periods to simulate
	double dt; // Time step size
	int t_points; // Number of time points to simulate
	double t_spacing; // Time spacing for data output, isntead of logging each timestep
	double forcing_period; // The period of the environmental switching (forcing)

	//Space variables
	double L; // The box side length
	int N; // The number of grid points in each direction (N x N grid)
	double dx; // The spacing between grid points in real space

	//Systems Parameters
	double astar, omega0;
	double lambda(double t); // The predation rate depends on time
	double k(double t); // This is the forcing function which is 1 / carrying capacity (called \kappa in draft)
	double k0; // Average environment
	double k1; // Forcing amplitude
	double n; // Forcing period parameter (see draft for exact definition)
	double D; // Diffusion coefficient
	double alpha; // Homotopy parameter

	// Density arrays
	fftw_complex* a; // This is the predator density array
	fftw_complex* b; // This is the prey density array
	
	// ETD2 variables
	fftw_complex* f; // This is the non-linear term for the predator equation
	fftw_complex* g; // This is the non-linear term for the prey equation
    fftw_complex* f_prev; // This is the non-linear term of the predator equation at the previous time step
    fftw_complex* g_prev; // This is the non-linear term of the prey equation at the previous time step
	std::vector<double> c_a; // Linear coefficient for the predator equation 
	std::vector<double> c_b; // Linear coefficient for the prey equation
	std::vector<double> f1; // First ETD2 coefficient for the predator equation
	std::vector<double> g1; // First ETD2 coefficient for the prey equation
	std::vector<double> f2; // Second ETD2 coefficient for the predator equation
	std::vector<double> g2; // Second ETD2 coefficient for the prey equation
	std::vector<double> f3; // Third ETD2 coefficient for the predator equation
	std::vector<double> g3; // Third ETD2 coefficient for the prey equation
	std::vector<double> f4; // Fourth ETD2 coefficient for the predator equation
	std::vector<double> g4; // Fourth ETD2 coefficient for the prey equation

	// Space variables
	std::vector<double> xx; // x-space array
	std::vector<double> yy; // y-space array
	std::vector<double> kx; // kx-space array
	std::vector<double> ky; // ky-space array

	// FFTW specific variables used to peroform FFT
	fftw_plan forward_transform;
	fftw_plan backward_transform;
};

#endif

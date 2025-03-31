#ifndef ETD2_H
#define ETD2_H
#include <fftw3.h>
#include <string>
#include <vector>
#include <random>

class ETD2_solver {
	public:
		ETD2_solver(double t_start_init, int t_end_init, double L_init, int N_init, double ustar_init, double k0_init, double k1_init, double n_init, double alpha_init, double D_init, double tt_spacing);
		double complex_magnitude(fftw_complex z);
		void complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result);
		void solve_in_log();
		void solve_in_regular();
		~ETD2_solver();

	private:
	//Methods
	void calculate_wavenumbers();
	void calculate_coeffcients_regular();
	void calculate_coeffcients_log();
	void set_initial_conditions(std::string type);
	void calculate_non_linear_log(double t);
	void calculate_non_linear_regular();
	void time_step(bool first_time);
	void write_data(std::string file_name, bool k_space = false);
	void transform_forward();
	void transform_to_k_space();
	void transform_to_x_space();
	void transform_to_log_space();
	void transform_to_regular_space();
	double calculate_omega0();

	//Variables
	//Time variables
	double t_start;
	double t_end;
	int num_of_periods;
	double dt;
	int t_points;
	double t_spacing;
	double forcing_period;

	//Space variables
	double L;
	int N;
	double dx;

	//Systems Parameters
	double ustar, omega0;
	double lambda(double t);
	double k(double t);
	double k0;
	double k1;
	double n;
	double D;
	double alpha;

	fftw_complex* u;
	fftw_complex* v;
	fftw_complex* f;
	fftw_complex* g;
    fftw_complex* f_prev;
    fftw_complex* g_prev;
	std::vector<double> c_u;
	std::vector<double> c_v;
	std::vector<double> f1;
	std::vector<double> g1;
	std::vector<double> f2;
	std::vector<double> g2;
	std::vector<double> f3;
	std::vector<double> g3;
	std::vector<double> f4;
	std::vector<double> g4;
	std::vector<double> xx;
	std::vector<double> yy;
	std::vector<double> kx; 
	std::vector<double> ky; 
	fftw_plan forward_transform;
	fftw_plan backward_transform;
};

#endif

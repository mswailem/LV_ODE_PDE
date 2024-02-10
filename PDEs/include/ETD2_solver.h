#ifndef ETD2_H
#define ETD2_H
#include <fftw3.h>
#include <string>
#include <vector>

class ETD2_solver {
	public:
		ETD2_solver(double t_start_init, double t_end_init, double dt_init, double L_init, int N_init, double mu_init, double sigma_init, double lambda_init, double k_init, double Dx_init, double Dy_init);
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
	void set_initial_conditions(int type);
	void calculate_non_linear_log();
	void calculate_non_linear_regular();
	void time_step(bool first_time);
	void write_data(std::string file_name, bool k_space = false);
	void transform_forward();
	void transform_to_k_space();
	void transform_to_x_space();
	void transform_to_log_space();
	void transform_to_regular_space();

	//Variables
	//Time variables
	double t_start;
	double t_end;
	double dt;
	int t_points;

	//Space variables
	double L;
	int N;
	double dx;

	//Systems Parameters
	double mu;
	double sigma;
	double lambda;
	double k;
	double Dx;
	double Dy;

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

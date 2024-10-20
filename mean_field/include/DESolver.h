#ifndef DESOLVER_H
#define DESOLVER_H
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

// For debugging, otherwise comment out
#include <iostream>

struct Params {
	double us, k0, k1, alpha, n, wn, du, dv;
};

class DESolver
{
	public:

		DESolver(int Type);
		~DESolver();
		void solve(double t_initial, double t_final);
		std::vector<double> get_y();
		double get_period();
		void set_y(double Y0, double Y1) { y[0] = Y0; y[1] = Y1; };
		void set_k0(double K0) { params.k0 = K0; }
		void set_k1(double K1) { params.k1 = K1; }
		void set_wavenumber(double k) { params.wn = k; }
		void set_diff(double d_u, double d_v) {params.du = d_u; params.dv = d_v;}
		void set_n(double N) { params.n = N; }
		void set_alpha(double Alpha) { params.alpha = Alpha; }
		void set_fp(double Ustar) { params.us = Ustar;}
		void set_params(std::unordered_map<std::string, double> p);
		void initialize();
		Params get_params() { return params; }
		bool check_param_region();

	private:
		double calculate_omega0();
		double lambda(double k);
		double k(double t);
		static int full_eq (double t, const double y[], double f[], void *params); //This function is static because gsl requires it to be
		static int linear_eq (double t, const double y[], double f[], void *params); //This function is static because gsl requires it to be
		double get_max_k0();
		double get_max_k1();

		int type; //This is the type of equation to solve. 0 for full, 1 for linear
		Params params;
		double omega0, lambda0;
		double y[2];
		gsl_odeiv2_system sys;
		gsl_odeiv2_driver *d;

	// Functions for debugging
	public:
	void d_print_omega0() { std::cout << omega0 << std::endl; }
	void d_print_lambda0() { std::cout << lambda0 << std::endl; }
	void d_print_params() {std::cout << "us: " << params.us << " k0: " << params.k0 << " k1: " << params.k1 << " alpha: " << params.alpha << " n: " << params.n << " wn: " << params.wn << " du: " << params.du << " dv: " << params.dv << std::endl;}
	void d_print_period() {std::cout << get_period() << std::endl;}
	void d_print_k(double t) {std::cout << k(t) << std::endl;}
	void d_print_lambda(double t) {std::cout << lambda(k(t)) << std::endl;}

};
#endif //DESOLVER_H

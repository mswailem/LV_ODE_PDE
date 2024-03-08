#ifndef DESOLVER_H
#define DESOLVER_H
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>

class DESolver
{
	public:
		DESolver(double Ustar, double Vstar, double K0, double K1, double N, double Alpha, double Y0, double Y1, int Type);
		DESolver(int Type);
		~DESolver();
		void solve(double t_initial, double t_final, double dt);
		std::vector<double> get_y();
		double get_period();
		void set_y(double Y0, double Y1); 
		void set_k0(double K0) { k0 = K0; }
		void set_k1(double K1) { k1 = K1; }
		void set_n(double N) { n = N; }
		void set_alpha(double Alpha) { alpha = Alpha; }
		void set_ustar(double Ustar) { ustar = Ustar; }
		void set_vstar(double Vstar) { vstar = Vstar; }
		void initialize(); //This is used if the constructor without arguments is used

	private:
		double calculate_omega0();
		double lambda();
		double k();
		double mu();
		static int full_eq (double t, const double y[], double f[], void *params); //This function is static because gsl requires it to be
		static int linear_eq (double t, const double y[], double f[], void *params); //This function is static because gsl requires it to be

		int type; //This is the type of equation to solve. 0 for full, 1 for linear
		double t, ustar, vstar, k0, k1, n, alpha;
		double omega0, lambda0, mu0;
		double y[2];
		gsl_odeiv2_system sys;
		gsl_odeiv2_driver *d;

};
#endif //DESOLVER_H

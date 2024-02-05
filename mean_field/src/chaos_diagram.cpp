#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

//Calculating the intrinsic frequency
double omega0(double a0, double b0, double k) {
	double term1 = (b0/a0)*(b0*k-1)*(b0*k-1);
	double term2 = 0.25*b0*k*(-4+3*b0*k);
	return sqrt(term1+term2);
}

double lambda(double a0, double b0, double k, double k0,double alpha) {
	double lambda0 = (1/a0)*(1-(b0*k0))-(k0);
	double lambda1 = (1/a0)*(1-(b0*k))-(1*k);
	return (1-alpha)*lambda0+alpha*lambda1;
}

double k(double t, double k0, double k1, double n, double a0, double b0) {
	return k0+k1*cos((1/n)*t*omega0(a0,b0,k0));
}

double mu(double a0,double b0, double k, double k0, double alpha) {
	double mu0 = (b0/a0)*(1-(b0*k0))-(b0*k0);
	double mu1 = (b0/a0)*(1-(b0*k))-(b0*k);
	return (1-alpha)*mu0+alpha*mu1;
}

double sigma() {
	return 1;
}

//Differential equation
int func (double t, const double y[], double f[], void *params) {
	double* ps = (double *)params;
	double a0 = ps[0];
	double b0 = ps[1];
	double k0 = ps[2];
	double k1 = ps[3];
	double n = ps[4];
	double alpha = ps[5];
	f[0] = y[0]*(lambda(a0,b0,k(t,k0,k1,n,a0,b0),k0,alpha)*y[1]-mu(a0,b0,k(t,k0,k1,n,a0,b0),k0,alpha));
	f[1] = y[1]*(sigma()-(y[0]+y[1])*k(t,k0,k1,n,a0,b0)-lambda(a0,b0,k(t,k0,k1,n,a0,b0),k0,alpha)*y[0]);
    return GSL_SUCCESS;
}


int main(int argc, char* argv[]) {

	//Checking command line arugments
	if (argc < 11) {
		std::cout << "Required arguments: t_start points a0 b0 k0 k1 n alpha_spacing initial_conditions" << std::endl;
		return 0;
	}

	double* ps = new double[6];
    const double t_start = std::stod(argv[1]);
    const int points = std::atoi(argv[2]);
	ps[0] = std::stod(argv[3]); // a0
	ps[1] = std::stod(argv[4]); // b0
	ps[2] = std::stod(argv[5]); // k0
	ps[3] = std::stod(argv[6]); // k1
	ps[4] = std::stod(argv[7]); // n
	double alpha_spacing = std::stod(argv[8]);

	/* if (K1 > K0 || K0 > 1/(2*(a0+b0))) { cout << "1 must be less than k0 and k0 has to be less than 1/(2*(a0+b0)) otherwise k becomes negative" << endl; return 0; } */

	double t_end = 2*M_PI*ps[4]/(omega0(ps[0],ps[1],ps[2]));
	double dt = t_end/points;

	char* filename; 
	filename = new char[100];
	sprintf(filename,"../output/chaos_diagram/a0=%g_b0=%g_k0=%g_k1=%g_n=%g.dat", ps[0], ps[1], ps[2], ps[3], ps[4]);
	std::ofstream outfile(filename);

	ps[5] = 1;
	while (ps[5] > 0) {
		std::cout << "alpha=" << ps[5] << std::endl;
		double t = 0;
		double y[2] = { std::stod(argv[9]), std::stod(argv[10]) };
		gsl_odeiv2_system sys = {func, NULL, 2, ps};
		gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

		double a_start = 0;
		double b_start = 0;
		double t_0 = 0;
		bool first_time = true;
		int counter = points;
		for (int i = 0; i < 50*points + std::floor(t_start/dt); i++) {
			if (t>=t_start) {
				if (first_time) {
					first_time = false;
					a_start = y[0];
					b_start = y[1];
					t_0 = i;
					outfile << ps[5] << " " << y[0] << std::endl;
				}
				if ( (i - t_0) == counter ) {
					if (abs(y[0]-a_start) < 1e-6 && abs(y[1]-b_start) < 1e-6) {
						break;
					}
					outfile << ps[5] << " " << y[0] << std::endl;
					counter+=points;
				}
			}
			int status = gsl_odeiv2_driver_apply (d, &t, t+dt, y);
			if (status != GSL_SUCCESS) {
				printf ("error, return value=%d\n", status);
				break;
			}
		}
		gsl_odeiv2_driver_free (d);
		ps[5] = ps[5]-alpha_spacing; // alpha
	}
	delete[] ps;
}

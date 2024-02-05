#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

//Calculates the intrinsic frequency
double omega0(double a0, double b0, double k) {
	double term1 = (b0/a0)*(b0*k-1)*(b0*k-1);
	double term2 = 0.25*b0*k*(-4+3*b0*k);
	return sqrt(term1+term2);
}

double lambda(double a0, double b0, double k) {
	return (1/a0)*(1-b0*k)-k;
}

double k(double t, double k0, double k1, double n, double a0, double b0) {
	return k0+k1*cos((1/n)*t*omega0(a0,b0,k0));
}

double mu(double b0, double lambda) {
	return b0*lambda;
}

double sigma() {
	return 1;
}

//Differential equations
int func (double t, const double y[], double f[], void *params) {
	double* ps = (double *)params;
	double a0 = ps[0];
	double b0 = ps[1];
	double k0 = ps[2];
	double k1 = ps[3];
	double n = ps[4];
	f[0] = y[0]*(lambda(a0,b0,k(t,k0,k1,n,a0,b0))*y[1]-mu(b0,lambda(a0,b0,k(t,k0,k1,n,a0,b0))));
	f[1] = y[1]*(sigma()-(y[0]+y[1])*k(t,k0,k1,n,a0,b0)-lambda(a0,b0,k(t,k0,k1,n,a0,b0))*y[0]);
    return GSL_SUCCESS;
}


int main(int argc, char* argv[]) {

	if (argc < 12) {
		std::cout << "Required arguments: t_start t_end delta_t a0 b0 k0 k1 n initial_conditions filename" << std::endl;
		return 0;
	}

	double* ps = new double[5];
    double t_start = std::stod(argv[1]);
    double t_end = std::stod(argv[2]);
    double dt = std::stod(argv[3]);
	ps[0] = std::stod(argv[4]); // a0
	ps[1] = std::stod(argv[5]); // b0
	ps[2] = std::stod(argv[6]); // k0
	ps[3] = std::stod(argv[7]); // k1
	ps[4] = std::stod(argv[8]); // n
	std::string filename = argv[11];
	/* if (K1 > K0 || K0 > 1/(2*(a0+b0))) { cout << "1 must be less than k0 and k0 has to be less than 1/(2*(a0+b0)) otherwise k becomes negative" << endl; return 0; } */

    const int points = (int)ceil((t_end-t_start)/dt);

	std::ofstream outfile("../output/solve/" + filename);

	double t = 0;
	double y[2] = { std::stod(argv[9]), std::stod(argv[10]) };
	gsl_odeiv2_system sys = {func, NULL, 2, ps};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

	for (int i = 0; i < points; i++) {
		outfile << t << " " << y[0] << " " << y[1] << "\n";
		int status = gsl_odeiv2_driver_apply (d, &t, t+dt, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
	}
	gsl_odeiv2_driver_free (d);
	delete[] ps;
}

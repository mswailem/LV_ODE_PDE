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
	double lambda1 = (1/a0)*(1-(b0*k))-k;
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

//Differential equations
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
	if (argc < 12) {
		std::cout << "Required arguments: t_start t_end delta_t a0 b0 k0 k1 n alpha initial_conditions" << std::endl;
		return 0;
	}
	double* ps = new double[6];
    double t_start = std::stod(argv[1]);
    double t_end = std::stod(argv[2]);
    double dt = std::stod(argv[3]);
	ps[0] = std::stod(argv[4]); // a0
	ps[1] = std::stod(argv[5]); // b0
	ps[2] = std::stod(argv[6]); // k0
	ps[3] = std::stod(argv[7]); // k1
	ps[4] = std::stod(argv[8]); // n
	ps[5] = std::stod(argv[9]); // alpha
	/* if (K1 > K0 || K0 > 1/(2*(a0+b0))) { cout << "1 must be less than k0 and k0 has to be less than 1/(2*(a0+b0)) otherwise k becomes negative" << endl; return 0; } */

    const int points = (int)ceil(t_end/dt);


	char* filename; 
	filename = new char[200];
	sprintf(filename,"output/a0=%g_b0=%g_k0=%g_k1=%g_n=%g_alpha=%g.dat", ps[0], ps[1], ps[2], ps[3], ps[4],ps[5]);
	std::ofstream outfile(filename);
	char* filename1; 
	filename1 = new char[200];
	sprintf(filename1,"output/a0=%g_b0=%g_k0=%g_k1=%g_n=%g_alpha=%g_fp.dat", ps[0], ps[1], ps[2], ps[3], ps[4],ps[5]);
	std::ofstream outfile1(filename1);

	double t = 0;
	double y[2] = { std::stod(argv[10]), std::stod(argv[11]) };
	gsl_odeiv2_system sys = {func, NULL, 2, ps};
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

	for (int i = 0; i < points; i++) {
		if (t > t_start) {
			outfile << t << " " << y[0] << " " << y[1] << "\n";
		}
		double kk = k(t,ps[2],ps[3],ps[4],ps[0],ps[1]);
		double lambdal = lambda(ps[0],ps[1],kk,ps[2],ps[5]);
		double mum = mu(ps[0],ps[1],kk,ps[2],ps[5]);
		double a = (1/(lambdal+kk))*(1-(mum*kk)/lambdal);
		double b = mum/lambdal;
		if (t > t_start) {
			outfile1 << t << " " << a << " " << b << std::endl;
		}
		int status = gsl_odeiv2_driver_apply (d, &t, t+dt, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
	}
	gsl_odeiv2_driver_free (d);
	delete ps;
}

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

//Calculating the intrinsic frequency
double omega0(double a0, double b0, double k) {
	double term1 = (b0/a0)*(b0*k-1)*(b0*k-1);
	double term2 = 0.25*b0*k*(-4+3*b0*k);
	return sqrt(term1+term2);
}

double k(double t, double k0, double k1, double n, double a0, double b0) {
	return k0+k1*cos((1/n)*t*omega0(a0,b0,k0));
}

//Linearezed ODEs
int func (double t, const double y[], double f[], void *params) {
	double* ps = (double *)params;
	double a0 = ps[0];
	double b0 = ps[1];
	double k0 = ps[2];
	double k1 = ps[3];
	double omega = ps[4];
	f[0] = -(((a0+b0)*(k(t,k0,k1,omega,a0,b0)))-1)*y[1];
	f[1] = -(b0/a0)*((a0*y[1]-b0*y[0])*k(t,k0,k1,omega,a0,b0)+y[0]);
	return GSL_SUCCESS;
}


int main(int argc, char* argv[]) {
	if (argc < 4) {
		std::cout << "Required arguments: points n b0_step" << std::endl;
		return 0;
	}
	double* ps = new double[5];
	double t_start = 0;
	int points = stoi(argv[1]);
	ps[4] = stod(argv[2]); // n
	const double b0_step = stod(argv[3]);
	/* if (K1 > K0 || K0 > 1/(2*(a0+b0))) { cout << "k1 must be less than k0 and k0 has to be less than 1/(2*(a0+b0)) otherwise k becomes negative" << endl; return 0; } */

	// Variable defintions
	double t_end;
	double dt;
	ofstream output("output/n="+string(argv[2])+".dat");

	//Allocating memory for gsl variables
	gsl_matrix *fundemental_matrix = gsl_matrix_alloc(2, 2);
	gsl_vector_complex *floquet_multipliers = gsl_vector_complex_alloc(2);
	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);

	for (int j = 0; j < 50; j++) {

			ps[0] = 0.1+j*0.1; // a0
			ps[1] = b0_step; // b0
							 
			while (ps[1] <= ps[0]) {

				ps[2] = 0.98/(2*(ps[1]+ps[0])); // k0
				ps[3] = 0.98*ps[2]; // k1
				t_end = 2*M_PI*ps[4]/(omega0(ps[0],ps[1],ps[2]));
				dt = t_end/points;

				cout << "a0: " << ps[0] << " b0: " << ps[1] << endl;
				
				//Initial conditions
				double t = 0;
				double y[2] = {1,0};

				gsl_odeiv2_system sys = {func, NULL, 2, ps};
				gsl_odeiv2_driver *diff = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

				//Solving the system
				for (int i = 0; i < points; i++) {
					int status = gsl_odeiv2_driver_apply (diff, &t, t+dt, y);
					if (status != GSL_SUCCESS) {
						printf ("error, return value=%d\n", status);
						break;
					}
				}

				//Assigning the values of the fundemental matrix
				gsl_matrix_set(fundemental_matrix, 0, 0, y[0]);
				gsl_matrix_set(fundemental_matrix, 0, 1, y[1]);

				//Resetting initial conditions
				y[0] = 0;
				y[1] = 1;
				t = 0;

				//Solving the system
				for (int i = 0; i < points; i++) {
					int status = gsl_odeiv2_driver_apply (diff, &t, t+dt, y);
					if (status != GSL_SUCCESS) {
						printf ("error, return value=%d\n", status);
						break;
					}
				}

				//Assigning the values of the fundemental matrix
				gsl_matrix_set(fundemental_matrix, 1, 0, y[0]);
				gsl_matrix_set(fundemental_matrix, 1, 1, y[1]);

				//Computing the floquet multipliers
				gsl_eigen_nonsymm(fundemental_matrix, floquet_multipliers, w);

				if (gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,0)) > 1 || gsl_complex_abs(gsl_vector_complex_get(floquet_multipliers,1)) > 1) {
					output << ps[0] << " " << ps[1] << " " << ps[2] << " " << ps[3] << std::endl;
				}
				gsl_odeiv2_driver_free (diff);
				ps[1] += b0_step;
			}
	}

	//Freeing the memory
	delete ps;
	gsl_matrix_free(fundemental_matrix);
	gsl_vector_complex_free(floquet_multipliers);
	gsl_eigen_nonsymm_free(w);
}

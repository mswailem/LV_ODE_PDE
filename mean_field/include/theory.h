#ifndef THEORY_H
#define THEORY_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex_double.h>

inline gsl_vector_complex* get_eigenvalues(double d_u, double d_v, double u0, double v0, double k0, double wn) {

	gsl_matrix *Jacobian = gsl_matrix_alloc(2, 2);
	gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(2);
	
	double a = -d_u*wn*wn;
	double b = -(u0+v0)*k0+1;
	double c = (v0*v0/u0)*k0-(v0/u0);
	double d = -d_v*wn*wn-v0*k0;

	gsl_matrix_set(Jacobian, 0, 0, a);
	gsl_matrix_set(Jacobian, 0, 1, b);
	gsl_matrix_set(Jacobian, 1, 0, c);
	gsl_matrix_set(Jacobian, 1, 1, d);

	gsl_vector_complex *eigenvalues = gsl_vector_complex_alloc(2);
	
	gsl_eigen_nonsymm(Jacobian, eigenvalues, w);

	gsl_matrix_free(Jacobian);
	gsl_eigen_nonsymm_free(w);
	return eigenvalues;
}


#endif //THEORY_H

#include "DESolver.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <cmath>
#include <iostream> //For debugging


DESolver::DESolver(double Ustar, double Vstar, double K0, double K1, double N, double Alpha, double Y0, double Y1, int Type) :
	ustar(Ustar), vstar(Vstar), k0(K0), k1(K1), n(N), alpha(Alpha), type(Type)
{
	omega0 = calculate_omega0();
	lambda0 = (1/ustar)*(1-(vstar*k0))-(k0); 
	mu0 = vstar * lambda0;
	y[0] = Y0;
	y[1] = Y1;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = this;
	if (type == 0) {
		sys.function = full_eq;
	} else {
		sys.function = linear_eq;
	}
	d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
}

void DESolver::initialize() {
	omega0 = calculate_omega0();
	lambda0 = (1/ustar)*(1-(vstar*k0))-(k0); 
	mu0 = vstar * lambda0;
}

DESolver::DESolver(int Type) :
	type(Type)
{
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = this;
	if (type == 0) {
		sys.function = full_eq;
	} else {
		sys.function = linear_eq;
	}
	d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
}

DESolver::~DESolver() {
	gsl_odeiv2_driver_free(d);
}

std::vector<double> DESolver::get_y() { return {y[0], y[1]}; }

void DESolver::set_y(double Y0, double Y1) { y[0] = Y0; y[1] = Y1; }

double DESolver::get_period() { return 2*M_PI/((1/n)*omega0); }

void DESolver::solve(double t_initial, double t_final, double dt) {
	t = t_initial;
	while (t < t_final) {
		int status = gsl_odeiv2_driver_apply(d, &t, t+dt, y);
		if (status != GSL_SUCCESS) {
			break;
		}
	}
}

double DESolver::calculate_omega0() {
	double term1 = (vstar/ustar)*(vstar*k0-1)*(vstar*k0-1);
	double term2 = 0.25*vstar*k0*(-4+3*vstar*k0);
	return sqrt(term1+term2);
}

double DESolver::lambda() {
	double lambda1 = (1/ustar)*(1-(vstar*k()))-k();
	return (1-alpha)*lambda0+alpha*lambda1;
}

double DESolver::k() {
	if (n == 0) {
		return k0;
	}
	return k0+k1*cos((1/n)*t*omega0);
}

double DESolver::mu() {
	double mu1 = (vstar/ustar)*(1-(vstar*k()))-(vstar*k());
	return (1-alpha)*mu0+alpha*mu1;
}

int DESolver::full_eq(double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	f[0] = y[0]*(solver->lambda()*y[1]-solver->mu());
	f[1] = y[1]*(1-(y[0]+y[1])*solver->k()-solver->lambda()*y[0]);
    return GSL_SUCCESS;
}

int DESolver::linear_eq (double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	f[0] = -(((solver->ustar+solver->vstar)*(solver->k()))-1)*y[1];
	f[1] = -(solver->vstar/solver->ustar)*((solver->ustar*y[1]-solver->vstar*y[0])*solver->k()+y[0]);
	return GSL_SUCCESS;
}

#include "DESolver.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <cmath>
#include <iostream> //For debugging


// TODO: I am currently playing around with optimization of the error tolerance and timestepping of my ODE solver

void DESolver::initialize() {
	omega0 = calculate_omega0();
	lambda0 = (1/params.us)*(1-(params.vs*params.k0))-(params.k0); 
	mu0 = params.vs * lambda0;
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
	d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-2, 1e-4, 0.0);
}

DESolver::~DESolver() {
	gsl_odeiv2_driver_free(d);
}

void DESolver::set_params(std::unordered_map<std::string, double> p) {
	params.us = p["us"];
	params.vs = p["vs"];
	params.k0 = p["k0"];
	params.k1 = p["k1"];
	params.alpha = p["alpha"];
	params.n = p["n"];
	params.wn = p["wn"];
	params.du = p["du"];
	params.dv = p["dv"];
}

std::vector<double> DESolver::get_y() { return {y[0], y[1]}; }

double DESolver::get_period() { return 2*M_PI/((1/params.n)*omega0); }

void DESolver::solve(double t_initial, double t_final) {
	double t = t_initial;
	int status = gsl_odeiv2_driver_apply(d, &t, t_final, y);
	if (status != GSL_SUCCESS) {
		std::cerr << "Error: " << gsl_strerror(status) << std::endl;
	}
}

double DESolver::calculate_omega0() {
	double term1 = (params.vs/params.us)*(params.vs*params.k0-1)*(params.vs*params.k0-1);
	double term2 = 0.25*params.vs*params.k0*(-4+3*params.vs*params.k0);
	return sqrt(term1+term2);
}

double DESolver::lambda(double k) {
	double lambda1 = (1/params.us)*(1-(params.vs*k))-k;
	return (1-params.alpha)*lambda0+params.alpha*lambda1;
}

double DESolver::k(double t) {
	if (params.n == 0) {
		return params.k0;
	}
	return params.k0+params.k1*cos((1/params.n)*t*omega0);
}

double DESolver::mu(double l) {
	return params.vs*l;
}

int DESolver::full_eq(double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	double cc = solver->k(t);
	double l = solver->lambda(cc);
	double m = solver->mu(l);
	f[0] = y[0]*(l*y[1]-m);
	f[1] = y[1]*(1-(y[0]+y[1])*cc-l*y[0]);
	return GSL_SUCCESS;
}

int DESolver::linear_eq (double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	double cc = solver->k(t);
	double us = solver->params.us;
	double vs = solver->params.vs;
	double wn = solver->params.wn;
	double ddu = solver->params.du;
	double ddv = solver->params.dv;
	f[0] = -ddu*wn*wn*y[0]-(((us+vs)*cc)-1)*y[1]; //Implement diffusion constant
	f[1] = -ddv*wn*wn*y[1]-(vs/us)*((us*y[1]-vs*y[0])*cc+y[0]); //Implement diffusion constant
	return GSL_SUCCESS;
}

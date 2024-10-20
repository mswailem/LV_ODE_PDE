#include "DESolver.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <cmath>
#include <iostream> //For debugging

void DESolver::initialize() {
	omega0 = calculate_omega0();
	lambda0 = (1/params.us)*(1-params.k0) - params.k0;
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
	d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-10, 1e-10);
}

DESolver::~DESolver() {
	gsl_odeiv2_driver_free(d);
}

void DESolver::set_params(std::unordered_map<std::string, double> p) {
	params.us = p["us"];
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
	double term1 = 4*pow(1-params.k0, 2)-4*params.us*params.k0*(1-params.k0)-params.us*pow(params.k0,2);
	double term2 = 2*sqrt(params.us);
	return sqrt(term1)/term2;
}

double DESolver::lambda(double k) {
	double lambda1 = (1/params.us)*(1-k)-k;
	return (1-params.alpha)*lambda0+params.alpha*lambda1;
}

double DESolver::k(double t) {
	if (params.n == 0) {
		return params.k0;
	}
	return params.k0+params.k1*cos((1/params.n)*t*omega0);
}

int DESolver::full_eq(double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	double cc = solver->k(t);
	double l = solver->lambda(cc);
	f[0] = y[0]*l*(y[1]-1);
	f[1] = y[1]*(1-(y[0]+y[1])*cc-l*y[0]);
	// Handling round-off error, otherwise the solution at the fixed point might move away from the fixed point
	f[0] = abs(f[0]) < 1e-10 ? 0 : f[0];
	f[1] = abs(f[1]) < 1e-10 ? 0 : f[1];

	return GSL_SUCCESS;
}

int DESolver::linear_eq (double t, const double y[], double f[], void *params) {
	DESolver* solver = static_cast<DESolver*>(params); //This is a pointer to the object that called this function
	double us = solver->params.us;
	double wn = solver->params.wn;
	double ddu = solver->params.du;
	double ddv = solver->params.dv;
	double n = solver->params.n;
	double o0 = solver->omega0;
	double k0 = solver->params.k0;
	double k1 = solver->params.k1;
	double alpha = solver->params.alpha;
	double cos_term = cos((1/n)*t*o0);
	double a_term = 1-k0-k1*(1+(1+us)*(alpha-1))*cos_term;
	double b_term = k0+k1*(1-(1+us)*(alpha-1))*cos_term;
	double last_term = k1*(1+us)*(alpha-1)*cos_term;
	f[0] = -ddu*wn*wn*y[0]+y[1]*(1-(1+us)*(k0+alpha*k1*cos_term));
	f[1] = -ddv*wn*wn*y[1]-(y[0]/us)*a_term-y[1]*b_term+last_term;
	return GSL_SUCCESS;
}

double DESolver::get_max_k0() {
	double term1 = params.us*(1+params.us);
	double term2 = params.us+2-sqrt(term1);
	double term3 = 3*params.us+4;
	return 2*term2/term3;
}

double DESolver::get_max_k1() {
	double term1 = 1-(1+params.us)*params.k0;
	double term2 = params.alpha*(1+params.us);
	return std::min(term1/term2, params.k0);
}

bool DESolver::check_param_region() {
	if (params.k0 < 0 || params.k0 > get_max_k0()) {
		return false;
	}
	if (params.k1 < 0 || params.k1 > get_max_k1()) {
		return false;
	}
	if (params.alpha < 0 || params.alpha > 1) {
		return false;
	}
	if (params.n < 0) {
		return false;
	}
	if (params.us < 0) {
		return false;
	}
	return true;
}

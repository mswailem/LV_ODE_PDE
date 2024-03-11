#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include "ETD2_solver.h"

ETD2_solver::ETD2_solver(double t_start_init, double t_end_init, double dt_init, double L_init, int N_init, double ustar_init, double vstar_init, double k0_init, double k1_init, double n_init, double Dx_u_init, double Dy_u_init, double Dx_v_init, double Dy_v_init, int tt_spacing) :
	t_start(t_start_init), t_end(t_end_init), dt(dt_init), L(L_init), N(N_init), ustar(ustar_init), vstar(vstar_init),  k0(k0_init), k1(k1_init), n(n_init), Dx_u(Dx_u_init), Dy_u(Dy_u_init), Dx_v(Dx_v_init), Dy_v(Dy_v_init), t_spacing(tt_spacing),
	c_u(N*N), c_v(N*N), f1(N*N), g1(N*N), f2(N*N), g2(N*N), f3(N*N), g3(N*N), f4(N*N), g4(N*N), xx(N), yy(N), kx(N), ky(N)
{
	omega0 = calculate_omega0();
	dx = L / N;
	t_points = std::ceil((t_end - t_start) / dt);
    u = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    v = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    g = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    f_prev = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    g_prev = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);

	forward_transform = fftw_plan_dft_2d(N, N, u, u, FFTW_FORWARD, FFTW_MEASURE);
	backward_transform = fftw_plan_dft_2d(N, N, u, u, FFTW_BACKWARD, FFTW_MEASURE);

}

ETD2_solver::~ETD2_solver() {
	fftw_destroy_plan(forward_transform);
	fftw_destroy_plan(backward_transform);
	fftw_free(u);
	fftw_free(v);
	fftw_free(f);
	fftw_free(g);
	fftw_free(f_prev);
	fftw_free(g_prev);
}

double ETD2_solver::calculate_omega0() {
	double term1 = (vstar/ustar)*(vstar*k0-1)*(vstar*k0-1);
	double term2 = 0.25*vstar*k0*(-4+3*vstar*k0);
	return sqrt(term1+term2);
}

double ETD2_solver::lambda(double t) {
	return (1/ustar)*(1-(vstar*k(t)))-k(t);
}

double ETD2_solver::mu(double t) {
	return (vstar/ustar)*(1-(vstar*k(t)))-(vstar*k(t));
}

double ETD2_solver::k(double t) {
	if (n == 0) {
		return k0;
	}
	return k0+k1*cos((1/n)*t*omega0);
}

void ETD2_solver::calculate_wavenumbers() {
    double dk = 2 * M_PI / L;  // Wavenumber spacing
    int N_half = N / 2;

    for (int i = 0; i < N; i++) {
        // Adjust for negative frequencies
        kx[i] = (i < N_half) ? i * dk : (i - N) * dk;
        ky[i] = (i < N_half) ? i * dk : (i - N) * dk;
    }

}

void ETD2_solver::calculate_coeffcients_log() {
	double k2_u;
	double k2_v;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			k2_u = -Dx_u * kx[i] * kx[i] - Dy_u * ky[j] * ky[j];
			k2_v = -Dx_v * kx[i] * kx[i] - Dy_v * ky[j] * ky[j];
			c_u[i * N + j] = k2_u;
			f1[i * N + j] = exp(k2_u * dt);
			f2[i * N + j] = c_u[i * N + j] == 0 ? dt : (f1[i * N + j] - 1) / c_u[i * N + j];
			f3[i * N + j] = c_u[i * N + j] == 0 ? 1.5*dt : ((1+c_u[i * N + j]*dt)*f1[i * N + j] - 1 - 2*c_u[i * N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
			f4[i * N + j] = c_u[i * N + j] == 0 ? -0.5*dt : (-f1[i * N + j] + 1 + c_u[i* N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
			c_v[i * N + j] = k2_v;
			g1[i * N + j] = exp(k2_v * dt);
			g2[i * N + j] = c_v[i * N + j] == 0 ? dt : (g1[i * N + j] - 1) / c_v[i * N + j];
			g3[i * N + j] = c_v[i * N + j] == 0 ? 1.5*dt : ((1+c_v[i * N + j]*dt)*g1[i * N + j] - 1 - 2*c_v[i * N + j]*dt) / (std::pow(c_v[i * N + j],2)*dt);
			g4[i * N + j] = c_v[i * N + j] == 0 ? -0.5*dt : (-g1[i * N + j] + 1 + c_v[i* N + j]*dt) / (std::pow(c_v[i * N + j],2)*dt);
		}
	}
}

void ETD2_solver::calculate_coeffcients_regular() { //Need to adjust this for now mu depends on time
	double k2_u;
	double k2_v;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			k2_u = -Dx_u * kx[i] * kx[i] - Dy_u * ky[j] * ky[j];
			k2_v = -Dx_v * kx[i] * kx[i] - Dy_v * ky[j] * ky[j];
			c_u[i * N + j] = k2_u - mu(0);
			f1[i * N + j] = exp(c_u[i * N + j] * dt);
			f2[i * N + j] = c_u[i * N + j] == 0 ? dt : (f1[i * N + j] - 1) / c_u[i * N + j];
			f3[i * N + j] = c_u[i * N + j] == 0 ? 1.5*dt : ((1+c_u[i * N + j]*dt)*f1[i * N + j] - 1 - 2*c_u[i * N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
			f4[i * N + j] = c_u[i * N + j] == 0 ? -0.5*dt : (-f1[i * N + j] + 1 + c_u[i* N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
			c_v[i * N + j] = k2_v + 1;
			g1[i * N + j] = exp(c_v[i * N + j] * dt);
			g2[i * N + j] = c_v[i * N + j] == 0 ? dt : (g1[i * N + j] - 1) / c_v[i * N + j];
			g3[i * N + j] = c_v[i * N + j] == 0 ? 1.5*dt : ((1+c_v[i * N + j]*dt)*g1[i * N + j] - 1 - 2*c_v[i * N + j]*dt) / (std::pow(c_v[i * N + j],2)*dt);
			g4[i * N + j] = c_v[i * N + j] == 0 ? -0.5*dt : (-g1[i * N + j] + 1 + c_v[i* N + j]*dt) / (std::pow(c_v[i * N + j],2)*dt);
		}
	}
}

double ETD2_solver::complex_magnitude(fftw_complex z) {
	return sqrt(z[0] * z[0] + z[1] * z[1]);
}

void ETD2_solver::complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result) {
	fftw_complex result_temp;
	result_temp[0] = z[0]*w[0]-z[1]*w[1];
	result_temp[1] = z[0]*w[1]+z[1]*w[0];
	result[0] = result_temp[0];
	result[1] = result_temp[1];
}

void ETD2_solver::calculate_non_linear_log(double t) {
	fftw_complex* deriv_u_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex* deriv_v_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex* deriv_u_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	fftw_complex* deriv_v_y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	transform_to_k_space();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			deriv_u_x[i * N + j][0] =  kx[i] * u[i * N + j][1];
			deriv_u_x[i * N + j][1] = -kx[i] * u[i * N + j][0];
			deriv_v_x[i * N + j][0] =  kx[i] * v[i * N + j][1];
			deriv_v_x[i * N + j][1] = -kx[i] * v[i * N + j][0];
			deriv_u_y[i * N + j][0] =  ky[j] * u[i * N + j][1];
			deriv_u_y[i * N + j][1] = -ky[j] * u[i * N + j][0];
			deriv_v_y[i * N + j][0] =  ky[j] * v[i * N + j][1];
			deriv_v_y[i * N + j][1] = -ky[j] * v[i * N + j][0];
		}
	}
	transform_to_x_space();
	fftw_execute_dft(backward_transform, deriv_u_x, deriv_u_x);
	fftw_execute_dft(backward_transform, deriv_v_x, deriv_v_x);
	fftw_execute_dft(backward_transform, deriv_u_y, deriv_u_y);
	fftw_execute_dft(backward_transform, deriv_v_y, deriv_v_y);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f[i * N + j][0] = Dx_u*std::pow(deriv_u_x[i * N + j][0], 2) + Dy_u*std::pow(deriv_u_y[i * N + j][0], 2) +  lambda(t) * exp(v[i * N + j][0]) - mu(t);
			g[i * N + j][0] = Dx_v*std::pow(deriv_v_x[i * N + j][0], 2) + Dy_v*std::pow(deriv_v_y[i * N + j][0], 2) - lambda(t) * exp(u[i * N + j][0]) + (1-k(t)*exp(u[i * N +j][0])-k(t)*exp(v[i * N + j][0]));
			f[i * N + j][1] = 0;
			g[i * N + j][1] = 0;
		}
	}
	fftw_free(deriv_u_x);
	fftw_free(deriv_v_x);
	fftw_free(deriv_u_y);
	fftw_free(deriv_v_y);
}

void ETD2_solver::calculate_non_linear_regular() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f[i * N + j][0] = lambda(0) * v[i * N + j][0] * u[i * N + j][0];
			g[i * N + j][0] = -lambda(0) * v[i * N + j][0] * u[i * N + j][0] -(u[i * N + j][0] * v[i * N + j][0] + std::pow(v[i * N + j][0], 2));
			f[i * N + j][1] = 0;
			g[i * N + j][1] = 0;
		}
	}
}

void ETD2_solver::set_initial_conditions(int type) {
	for (int i = 0; i < N; i++) {
		xx[i] = i * dx;
		for (int j = 0; j < N; j++) {
			if (i == 0) {yy[j] = j * dx;}
			if (type == 0) { //Checkboard pattern (need to make type names more transparent)
				if ( xx[i] < (double)L / 2  && yy[j] < (double)L / 2) {
					u[i * N +j][0] = 1;
					u[i * N +j][1] = 0;
					v[i * N +j][0] = 0.01;
					v[i * N +j][1] = 0;
				} else if ( xx[i] < (double)L / 2 && yy[j] >= (double)L / 2) {
					u[i * N +j][0] = 0.01;
					u[i * N +j][1] = 0;
					v[i * N +j][0] = 1;
					v[i * N +j][1] = 0;
				} else if ( xx[i] >= (double)L / 2 && yy[j] < (double)L / 2) {
					u[i * N +j][0] = 0.01;
					u[i * N +j][1] = 0;
					v[i * N +j][0] = 1;
					v[i * N +j][1] = 0;
				} else {
					u[i * N +j][0] = 1;
					u[i * N +j][1] = 0;
					v[i * N +j][0] = 0.01;
					v[i * N +j][1] = 0;
				}
			} 
			else if (type == 1) { //Continuous initial conditions with sines and cosines (need to make type names more transparent)
				u[i * N + j][0] = exp(sin(2*M_PI*xx[i]/L));
				u[i * N + j][1] = 0;
				v[i * N + j][0] = exp(cos(2*M_PI*yy[j]/L));
				v[i * N + j][1] = 0;
			} else if (type == 2) { //Constant initial conditions (need to make type names more transparent)
				u[i * N + j][0] = ustar + 0.01;
				u[i * N + j][1] = 0;
				v[i * N + j][0] = vstar + 0.01;
				v[i * N + j][1] = 0;
			}

		}
	}
}

void ETD2_solver::write_data(std::string file_name, bool k_space) {
	if (k_space) {
		std::ofstream file("../output/" + file_name+"_k.dat");
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				file << kx[i] << " " << ky[j] << " " << complex_magnitude(u[i * N + j]) << " " << complex_magnitude(v[i * N +j]) << std::endl;
			}
			file << std::endl;
		}
		file.close();
	} else {
		std::ofstream file("../output/" + file_name+"_x.dat");
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				file << xx[i] << " " << yy[j] << " " << u[i * N + j][0] << " " << v[i * N +j][0] << std::endl;
			}
			file << std::endl;
		}
		file.close();
	}
}

void ETD2_solver::transform_to_k_space() {
	fftw_execute_dft(forward_transform, u, u);
	fftw_execute_dft(forward_transform, v, v);
	fftw_execute_dft(forward_transform, f, f);
	fftw_execute_dft(forward_transform, g, g);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			u[i * N + j][0] /= N*N;
			u[i * N + j][1] /= N*N;
			v[i * N + j][0] /= N*N;
			v[i * N + j][1] /= N*N;
			f[i * N + j][0] /= N*N;
			g[i * N + j][0] /= N*N;
			f[i * N + j][1] /= N*N;
			g[i * N + j][1] /= N*N;
		}
	}
}

void ETD2_solver::transform_to_x_space() {
	fftw_execute_dft(backward_transform, u, u);
	fftw_execute_dft(backward_transform, v, v);
	fftw_execute_dft(backward_transform, f, f);
	fftw_execute_dft(backward_transform, g, g);
}

void ETD2_solver::transform_to_log_space() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			u[i * N + j][0] = log(u[i * N + j][0]);
			u[i * N + j][1] = 0;
			v[i * N + j][0] = log(v[i * N + j][0]);
			v[i * N + j][1] = 0;
		}
	}
}

void ETD2_solver::transform_to_regular_space() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			u[i * N + j][0] = exp(u[i * N + j][0]);
			u[i * N + j][1] = 0;
			v[i * N + j][0] = exp(v[i * N + j][0]);
			v[i * N + j][1] = 0;
		}
	}
}

void ETD2_solver::time_step(bool first_time) {
	if (first_time) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				u[i * N + j][0] = f1[i * N + j] * u[i * N + j][0] + f2[i * N + j] * f[i * N + j][0];
				u[i * N + j][1] = f1[i * N + j] * u[i * N + j][1] + f2[i * N + j] * f[i * N + j][1];
				v[i * N + j][0] = g1[i * N + j] * v[i * N + j][0] + g2[i * N + j] * g[i * N + j][0];
				v[i * N + j][1] = g1[i * N + j] * v[i * N + j][1] + g2[i * N + j] * g[i * N + j][1];
				f_prev[i * N + j][0] = f[i * N + j][0];
				f_prev[i * N + j][1] = f[i * N + j][1];
				g_prev[i * N + j][0] = g[i * N + j][0];
				g_prev[i * N + j][1] = g[i * N + j][1];
			}
		}
	} else {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				u[i * N + j][0] = f1[i * N + j] * u[i * N + j][0] + f3[i * N + j] * f[i * N + j][0] + f4[i * N + j] * f_prev[i * N + j][0];
				u[i * N + j][1] = f1[i * N + j] * u[i * N + j][1] + f3[i * N + j] * f[i * N + j][1] + f4[i * N + j] * f_prev[i * N + j][1];
				v[i * N + j][0] = g1[i * N + j] * v[i * N + j][0] + g3[i * N + j] * g[i * N + j][0] + g4[i * N + j] * g_prev[i * N + j][0];
				v[i * N + j][1] = g1[i * N + j] * v[i * N + j][1] + g3[i * N + j] * g[i * N + j][1] + g4[i * N + j] * g_prev[i * N + j][1];
				f_prev[i * N + j][0] = f[i * N + j][0];
				f_prev[i * N + j][1] = f[i * N + j][1];
				g_prev[i * N + j][0] = g[i * N + j][0];
				g_prev[i * N + j][1] = g[i * N + j][1];
			}
		}
	}
}

void ETD2_solver::solve_in_log() {
	double t = t_start;
	calculate_wavenumbers();
	calculate_coeffcients_log();
	set_initial_conditions(1);
	write_data("0");
	transform_to_log_space();
	calculate_non_linear_log(t);
	transform_to_k_space();
	/* write_data("0", true); */
	int counter = t_spacing;
	for (int i = 0; i < t_points; i++) {
		t += dt;
		time_step(i == 0);
		transform_to_x_space();
		if (t >= counter) {
			std::cout << "Writing data at t = " << t << std::endl;
			transform_to_regular_space();
			write_data(std::to_string(counter));
			transform_to_log_space();
			counter += t_spacing;
		}
		calculate_non_linear_log(t);
		transform_to_k_space();
		/* if (t >= counter) { */
		/* 	write_data(std::to_string(counter), true); */
		/* 	counter++; */
		/* } */
	}
}

void ETD2_solver::solve_in_regular() {
	double t = t_start;
	calculate_wavenumbers();
	calculate_coeffcients_regular();
	set_initial_conditions(1);
	write_data("0");
	calculate_non_linear_regular();
	transform_to_k_space();
	/* write_data("0", true); */
	int counter = 1;
	for (int i = 0; i < t_points; i++) {
		t += dt;
		time_step(i == 0);
		transform_to_x_space();
		if (t >= counter) {
			write_data(std::to_string(counter));
			counter++;
		}
		calculate_non_linear_regular();
		transform_to_k_space();
		/* if (t >= counter) { */
		/* 	write_data(std::to_string(counter), true); */
		/* 	counter++; */
		/* } */
	}
}

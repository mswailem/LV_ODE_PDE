#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <fftw3.h>
#include "ETD2_solver.h"

ETD2_solver::ETD2_solver(double t_start_init, double t_end_init, double dt_init, double L_init, int N_init, double mu_init, double sigma_init, double lambda_init, double k_init, double Dx_init, double Dy_init) :
	t_start(t_start_init), t_end(t_end_init), dt(dt_init), L(L_init), N(N_init), mu(mu_init), sigma(sigma_init), lambda(lambda_init), k(k_init), Dx(Dx_init), Dy(Dy_init),
	c_u(N*N), c_v(N*N), f1(N*N), g1(N*N), f2(N*N), g2(N*N), f3(N*N), g3(N*N), f4(N*N), g4(N*N), xx(N), yy(N), kx(N), ky(N)
{
	
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

void ETD2_solver::calculate_wavenumbers() {
    double dk = 2 * M_PI / L;  // Wavenumber spacing
    int N_half = N / 2;

    for (int i = 0; i < N; i++) {
        // Adjust for negative frequencies
        kx[i] = (i < N_half) ? i * dk : (i - N) * dk;
        ky[i] = (i < N_half) ? i * dk : (i - N) * dk;
    }

}

void ETD2_solver::calculate_coeffcients() {
	double k2;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			k2 = -Dx * kx[i] * kx[i] - Dy * ky[j] * ky[j];
			c_u[i * N + j] = (k2 - mu);
			c_v[i * N + j] = (k2 + sigma);
			f1[i * N + j] = exp(c_u[i * N + j] * dt);
			f2[i * N + j] = c_u[i * N + j] == 0 ? dt : (f1[i * N + j] - 1) / c_u[i * N + j];
			f3[i * N + j] = c_u[i * N + j] == 0 ? 1.5*dt : ((1+c_u[i * N + j]*dt)*f1[i * N + j] - 1 - 2*c_u[i * N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
			f4[i * N + j] = c_u[i * N + j] == 0 ? -0.5*dt : (-f1[i * N + j] + 1 + c_u[i* N + j]*dt) / (std::pow(c_u[i * N + j],2)*dt);
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

void ETD2_solver::calculate_non_linear() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			f[i * N + j][0] = lambda * u[i * N + j][0] * v[i * N + j][0];
			g[i * N + j][0] = -sigma * (std::pow(v[i * N + j][0],2)/k) - lambda * u[i * N + j][0] * v[i * N + j][0];
			f[i * N + j][1] = 0;
			g[i * N + j][1] = 0;
		}
	}
}

void ETD2_solver::set_initial_conditions() {
	for (int i = 0; i < N; i++) {
		xx[i] = i * dx;
		for (int j = 0; j < N; j++) {
			if (i == 0) {yy[j] = j * dx;}
			if ( xx[i] < (double)L / 2  && yy[j] < (double)L / 2) {
				u[i * N +j][0] = 1;
				u[i * N +j][1] = 0;
				v[i * N +j][0] = 0;
				v[i * N +j][1] = 0;
			} else if ( xx[i] < (double)L / 2 && yy[j] >= (double)L / 2) {
				u[i * N +j][0] = 0;
				u[i * N +j][1] = 0;
				v[i * N +j][0] = 1;
				v[i * N +j][1] = 0;
			} else if ( xx[i] >= (double)L / 2 && yy[j] < (double)L / 2) {
				u[i * N +j][0] = 0;
				u[i * N +j][1] = 0;
				v[i * N +j][0] = 1;
				v[i * N +j][1] = 0;
			} else {
				u[i * N +j][0] = 1;
				u[i * N +j][1] = 0;
				v[i * N +j][0] = 0;
				v[i * N +j][1] = 0;
			}
			f[i * N + j][0] = lambda * u[i * N + j][0] * v[i * N + j][0];
			g[i * N + j][0] = -sigma * (std::pow(v[i * N + j][0],2)/k) - lambda * u[i * N + j][0] * v[i * N + j][0];
			f[i * N + j][1] = 0;
			g[i * N + j][1] = 0;
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
	} else {
		std::ofstream file("../output/" + file_name+"_x.dat");
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				file << xx[i] << " " << yy[j] << " " << u[i * N + j][0] << " " << v[i * N +j][0] << std::endl;
			}
			file << std::endl;
		}
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

void ETD2_solver::solve() {
	double t = t_start;
	calculate_wavenumbers();
	calculate_coeffcients();
	set_initial_conditions();
	write_data("0");
	transform_to_k_space();
	write_data("0", true);
	int counter = 1;
	for (int i = 0; i < t_points; i++) {
		t += dt;
		time_step(i == 0);
		transform_to_x_space();
		if (t >= counter) {
			write_data(std::to_string(counter));
			counter++;
		}
		calculate_non_linear();
		transform_to_k_space();
		if (t >= counter) {
			write_data(std::to_string(counter), true);
			counter++;
		}
	}
}

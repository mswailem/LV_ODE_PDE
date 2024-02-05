#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <fftw3.h>

//Function to calculate the wavenumbers for the FFT
std::vector<double> calculate_wavenumbers(int N, double L) {
    std::vector<double> wavenumbers(N);
    double dk = 2 * M_PI / L;  // Wavenumber spacing
    int N_half = N / 2;

    for (int i = 0; i < N; i++) {
        // Adjust for negative frequencies
        wavenumbers[i] = (i < N_half) ? i * dk : (i - N) * dk;
    }

    return wavenumbers;
}

//Function to precalculate the coefficients for the ETD2 method
void calculate_coeffcients(std::vector<double> kx, std::vector<double> ky, fftw_complex* c_u, fftw_complex* c_v, double mu, double sigma, double Dx, double Dy, double dt, int N) {
	double k2;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			k2 = -Dx * kx[i] * kx[i] - Dy * ky[j] * ky[j];
			c_u[i * N + j][0] = exp((k2 - mu) * dt);
			c_v[i * N + j][0] = exp((k2 + sigma) * dt);
			c_u[i * N + j][1] = 0;
			c_v[i * N + j][1] = 0;
		}
	}
}

//Function to calculate magnitude of complex number
double complex_magnitude(fftw_complex z) {
	return sqrt(z[0] * z[0] + z[1] * z[1]);
}

//Function that multiplies two complex numbers
void complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result) {
	fftw_complex result_temp;
	result_temp[0] = z[0]*w[0]-z[1]*w[1];
	result_temp[1] = z[0]*w[1]+z[1]*w[0];
	result[0] = result_temp[0];
	result[1] = result_temp[1];
}

//below is the nonlinear part of du/dt
double f(double u, double v,  double lambda) {
	return lambda * u * v;
}

//below is the nonlinear part of dv/dt
double g(double u, double v, double sigma, double k,  double lambda) {
	return -sigma * (std::pow(v,2)/k) - lambda * u * v;
}

int main(int argc, char** argv) {
	if (argc < 12) {
		std::cout << "Required arguments: t_start t_end dt L N mu sigma lambda kk Dx Dy" << std::endl;
		return 1;
	}
	
	//Reading arguments and defining variables
	double t_start = std::stod(argv[1]);
	double t_end = std::stod(argv[2]);
	double dt = std::stod(argv[3]);
	double L = std::stod(argv[4]);
	int N = std::stoi(argv[5]);
	double mu = std::stod(argv[6]);
	double sigma = std::stod(argv[7]);
	double lambda = std::stod(argv[8]);
	double k = std::stod(argv[9]);
	double Dx = std::stod(argv[10]);
	double Dy = std::stod(argv[11]);
	double dx = L / N;
	int t_points = std::ceil((t_end - t_start) / dt);

	//Defining arrays
    fftw_complex* u_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* v_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* u_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* v_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* f_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* g_x = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* f_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* g_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* c_u = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
    fftw_complex* c_v = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);
	
	//Initial conditions
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double x = i * dx;
			double y = j * dx;
			if ( x < L / 2  && y < L / 2) {
				u_x[i * N +j][0] = 1;
				u_x[i * N +j][1] = 0;
				v_x[i * N +j][0] = 0;
				v_x[i * N +j][1] = 0;
			} else if ( x < L / 2 && y >= L / 2) {
				u_x[i * N +j][0] = 0;
				u_x[i * N +j][1] = 0;
				v_x[i * N +j][0] = 1;
				v_x[i * N +j][1] = 0;
			} else if ( x >= L / 2 && y < L / 2) {
				u_x[i * N +j][0] = 0;
				u_x[i * N +j][1] = 0;
				v_x[i * N +j][0] = 1;
				v_x[i * N +j][1] = 0;
			} else {
				u_x[i * N +j][0] = 1;
				u_x[i * N +j][1] = 0;
				v_x[i * N +j][0] = 0;
				v_x[i * N +j][1] = 0;
			}
			f_x[i * N +j][0] = f(u_x[i * N +j][0], v_x[i * N +j][0], lambda);
			f_x[i * N +j][1] = 0;
			g_x[i * N +j][0] = g(u_x[i * N +j][0], v_x[i * N +j][0], sigma, k, lambda);
			g_x[i * N +j][1] = 0;
		}
	}


	//Below I need to precompute the coeffcients for the ETD2 method and then implement the ETDRK2 method

	return 0;
}


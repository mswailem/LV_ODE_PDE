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

//Function to calculate magnitude of complex number
double complex_magnitude(fftw_complex z) {
	return sqrt(z[0] * z[0] + z[1] * z[1]);
}

void complex_multiply(fftw_complex z, fftw_complex w, fftw_complex result) {
	fftw_complex result_temp;
	result_temp[0] = z[0]*w[0]-z[1]*w[1];
	result_temp[1] = z[0]*w[1]+z[1]*w[0];
	result[0] = result_temp[0];
	result[1] = result_temp[1];
}

int main() {
	int N = 64;
	double L = 10.0;
	double dx = L / N;
    fftw_complex* data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N);

	std::ofstream original_data("../output/data.dat");
	std::ofstream exact_derivative("../output/exact.dat");
	std::ofstream derivative_outfile("../output/derivative.dat");
	
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double x = i * dx;
			double y = j * dx;
			double value = sin(2 * M_PI * x / L) + cos(2 * M_PI * y / L);
			double derivative_value = -std::pow(2 * M_PI / L,2) * value; 
			data[i * N + j][0] = value;
			data[i * N + j][1] = 0;
			original_data << x << " " << y << " " << data[i * N +j][0] << std::endl;
			exact_derivative << x << " " << y << " " << derivative_value << std::endl;
		}
		original_data << std::endl;
		exact_derivative << std::endl;
	}

    // Perform 2D FFT
	fftw_plan forward_transform = fftw_plan_dft_2d(N, N, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan backward_transform = fftw_plan_dft_2d(N, N, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);
	std::vector<double> kx = calculate_wavenumbers(N, L);
	std::vector<double> ky = calculate_wavenumbers(N, L);
	fftw_execute(forward_transform);


	// Write FFT data to file
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fftw_complex k2;
			k2[0] =-kx[i] * kx[i] - ky[j] * ky[j];
			k2[1] = 0;
			data[i * N + j][0] *= 1.0/(N*N);
			data[i * N + j][1] *= 1.0/(N*N);
			complex_multiply(k2, data[i * N +j], data[i * N +j]);
		}
	}
	

	fftw_execute(backward_transform);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double x = i * dx;
			double y = j * dx;
			derivative_outfile << x << " " << y << " " << data[i * N + j][0] << std::endl;
		}
		derivative_outfile << std::endl;
	} 

	fftw_destroy_plan(forward_transform);
	fftw_destroy_plan(backward_transform);
	fftw_free(data);
	return 0;
}


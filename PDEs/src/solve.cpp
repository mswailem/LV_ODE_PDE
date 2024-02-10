#include <iostream>
#include <fftw3.h>
#include "ETD2_solver.h"


int main(int argc, char** argv) {
	
	if (argc < 12) {
		std::cout << "Required arguments: t_start t_end dt L N mu sigma lambda kk Dx Dy" << std::endl;
		return 1;
	}

	system("mkdir -p ../output");
	system("rm -f ../output/*");

	ETD2_solver solver(std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), std::stod(argv[10]), std::stod(argv[11]));

	solver.solve_in_regular();

	return 0;
}


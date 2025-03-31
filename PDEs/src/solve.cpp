#include "ETD2_solver.h"
#include <fftw3.h>
#include <iostream>

// TODO: Add a way of outputting the phase space plot for each pixel so that I can plot it
// TODO: Refactor the code to implement something similar to what I did with mean-field
int main(int argc, char **argv) {
  if (argc < 12) {
    std::cout << "Required arguments: t_start t_end L N ustar k0 k1 n alpha D t_spacing" << std::endl;
    return 1;
  }

  system("mkdir -p ../output");
  system("rm -f ../output/*");

  // Checking that we are in the physical regime
  double k0 = std::stod(argv[6]);
  double k1 = std::stod(argv[7]);
  double ustar = std::stod(argv[5]);
  if (k1 > k0 || k1 > 1 / (2 * (ustar + 1))) {
    std::cout << "Error: k0 must be smaller than k1 and k1 must be smaller "
                 "than 1/(2*(ustar + 1))"
              << std::endl;
    return 1;
  }

  ETD2_solver solver(std::stod(argv[1]), std::stoi(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stod(argv[9]), std::stod(argv[10]), std::stod(argv[11]));

  solver.solve_in_log();

  return 0;
}

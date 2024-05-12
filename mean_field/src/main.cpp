#include <iostream>
#include <vector>
#include <string>
#include "programs.h"
#include "iohandler.h"

// NOTE: Might change this so that instead each program has a function member variable
void run_program(Program program, std::unordered_map<std::string, double> params, std::vector<VaryingParam> vs, std::string filename) {
	if (program.name == "bifurcation diagram") {
		if (vs.size() == 1) {
			bifurcation_diagram(params, vs[0], params["t0"], filename);
		} else if (vs.size() == 2) {
			bifurcation_diagram(params, vs[0], vs[1], params["t0"], filename);
		} else {
			throw std::invalid_argument("Invalid number of varying parameters");
		}
	} else if (program.name == "time series") {
		time_series(params, params["t0"], params["tf"], params["dt"], {params["a0"], params["b0"]}, filename);
	} else if (program.name == "stability") {
		stability(params, vs[0], vs[1], filename);
	} else if (program.name == "dispersion relation") {
		dispersion_relation(params, vs[0], filename);
	} else {
		throw std::invalid_argument("Invalid program name");
	}
}

// NOTE: Need to implement a way to create params files and submit them maybe using a bash script? So make it optional if a param file is provided
int main(int argc, char* argv[]) {

	std::string picked_program_name = fzf_pick(get_program_names(), "Pick program: ");
	Program picked_program = get_program(picked_program_name);
	std::vector<VaryingParam> varying_params = pick_varying_params(picked_program);
	std::unordered_map<std::string, double> fixed_params = get_params(picked_program, varying_params);
	std::string filename = get_filename();
	run_program(picked_program, fixed_params, varying_params, filename);
	return 0;

}

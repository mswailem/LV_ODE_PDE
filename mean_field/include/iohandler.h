#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <iomanip>

// A struct to handle the logic of a varying parameter
struct VaryingParam {
	std::string name;
	double start;
	double end;
	double step;

	// Constructor
	VaryingParam(std::string name, double start, double end, double step) : name(name), start(start), end(end), step(step) {}

	// Get the number of points that this parameter will be varied over
	int get_num_of_points() {
		return ceil((end - start) / step);
	}
};

struct Program {
	std::string name;
	std::vector<std::string> required_params;
	int variable_params;
	int required_variable_params;
};

// NOTE: For testing, will implement the actual thing later
inline std::vector<Program> programs = {
	{"testing", {"t0", "points", "ustar", "k0", "n", "step"}, 2, 2},
	{"testing2", {"t0", "points", "ustar", "vstar", "k0", "n", "step"}, 2, 1}
};

// Get the names of all the programs
inline std::vector<std::string> get_program_names() {
	std::vector<std::string> names;
	for (auto const& x : programs) {
		names.push_back(x.name);
	}
	return names;
}

// Finds a program by name
inline Program get_program(std::string name) {
	for (auto const& x : programs) {
		if (x.name == name) {
			return x;
		} else {
			continue;
		}
	}
	throw std::invalid_argument("Invalid program name");
}

// Create the output file for a given program
inline std::ofstream create_outfile(std::string program, std::string filename) {
	std::string command = "mkdir -p ../output/" + program;
	std::system(command.c_str());
	std::ofstream out_file("../output/" + program + "/" + filename);
	return out_file;
}

// Formatting the double values to strings and ignoring trailiing zeros
inline std::string formatDouble(double value) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(6) << value; // Use fixed-point notation
    std::string str = stream.str();
    
    // Remove trailing zeros
    str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    
    // If the string ends with a decimal point, remove it as well
    if (str.back() == '.') {
        str.pop_back();
    }
    
    return str;
}

// Default filename based on the parameters of the program
inline std::string default_filename(std::unordered_map<std::string, double> p) {
    std::string filename = "";
    for (auto const& x : p) {
        filename += x.first + "=" + formatDouble(x.second) + "_";
    }
    filename.pop_back(); // Remove the last underscore
	filename += ".dat";
    return filename;
}

// Execute a command and return the output as a string
inline std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// Get user to pick from options using fzf
inline std::string fzf_pick(std::vector<std::string> options, const std::string prompt) {
	std::string command = "echo '";
	for (auto const& x : options) {
		command += x + "\n";
	}
	command.pop_back();
	command += "' | fzf --prompt='" + prompt + " > '";
	std::string choice = exec(command.c_str());
	choice.pop_back(); // Remove the newline character
	
	// Ensure that choice picked is in the options
	for (auto const& x : options) {
		if (choice == x) {
			return choice;
		} else {
			continue;
		}
	}
	throw std::invalid_argument("Invalid choice");
	return choice;
}

// Get the user to pick the parameters to vary
inline std::vector<VaryingParam> pick_varying_params(Program p) {
	std::vector<VaryingParam> varying_params;
	std::vector<std::string> options = p.required_params;
	std::string prompt = "Pick varying parameters ";
	for (int i = 0; i < p.variable_params; i++) {
		prompt += "v" + std::to_string(i+1);
		if (i < p.required_variable_params) {
			prompt += " (Required)";
		}
		prompt += " ";
	}
	prompt += ": \n Available options are: ";
	for (auto const& x : options) {
		prompt += x + " ";
	}
	std::cout << prompt << std::endl;
	 // Read a line of input from the user
    std::string inputLine;
	std::vector<std::string> choices;
    std::getline(std::cin, inputLine);
    
    std::istringstream iss(inputLine);
    std::string token;
    while (iss >> token) {
        choices.push_back(token);
    }
	for (auto const& x : choices) {
		double start, end, step;
		std::cout << "Enter the start, end and step (separated by spaces) for " << x << ": " << std::endl;
		std::cin >> start >> end >> step;
		varying_params.push_back(VaryingParam(x, start, end, step));
	}
	return varying_params;
}

#endif // IOHANDLER_H

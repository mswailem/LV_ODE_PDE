#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <iomanip>

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

#endif // IOHANDLER_H

#ifndef IOHANDLER_H
#define IOHANDLER_H

#include <iostream>
#include <fstream>
#include <string>

// TODO: Create a default naming convention based on the parameters of the program

// Create the output file for a given program
inline std::ofstream create_outfile(std::string program, std::string filename) {
	std::string command = "mkdir -p ../output/" + program;
	std::system(command.c_str());
	std::ofstream out_file("../output/" + program + "/" + filename);
	return out_file;
}

#endif // IOHANDLER_H

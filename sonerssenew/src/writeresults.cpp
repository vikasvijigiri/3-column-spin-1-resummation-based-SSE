#include <iostream>
#include <fstream>
#include <iomanip>


#include "../include/SSEvariables.hpp"
#include "../include/writeresults.hpp"

using namespace std;


void writeresults::save_bin_data(const char* filename , const char** headers , const double* values, int no_of_observables) {
    std::ofstream outputFile(filename, std::ios::app);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return;  // Return an error code or throw an exception if needed
    }

    // Check if the file is empty (to determine whether to write the header)
    outputFile.seekp(0, std::ios::end);
    bool fileIsEmpty = outputFile.tellp() == 0;

    // Write the header if the file is empty
    if (fileIsEmpty) {
        //const char* headers = createHeadersArray();

        // Write the header to the file
        for (int i = 0; i < no_of_observables; i++) {
                outputFile << std::setw(15) << headers[i];
        }

        outputFile << '\n' << '\n';
    }
    for (int i = 0; i < no_of_observables; i++) {
        outputFile << std::setw(15) << std::scientific << values[i];
    }


    outputFile << '\n';

    // Close the file
    outputFile.close();
}


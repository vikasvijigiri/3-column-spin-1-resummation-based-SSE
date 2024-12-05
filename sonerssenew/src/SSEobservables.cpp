#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>

#include "../include/SSEvariables.hpp"
#include "../include/SSEobservables.hpp"


struct Site {
    int index;
    int x, y;
};


// Function that initializes all the observables to zero (for measurement). 
void SSEobservables::measure_observables(SSEvariables& vars) {
    
    double E = -vars.n1/vars.Beta/vars.Ns/vars.iter;
    // Energy_density, Energy_density^2, Energy_density^4    
    enrg  += E;
    enrg2 += E*E;
    enrg4 += E*E*E*E;

    // Stiffness
    stiffx = 0.;
    stiffy = 0.;

}




// Function that initializes all the observables to zero (for measurement). 
void SSEobservables::Initiate_observables() {
    
    // Energy_density, Energy_density^2, Energy_density^4    
    enrg = 0.;
    enrg2 = 0.;
    enrg4 = 0.;

    // Stiffness
    stiffx = 0.;
    stiffy = 0.;

}


// 
const char** SSEobservables::createHeadersArray() {
    static const char* headers[5];

    // energy headers
    headers[0] = "enrg";
    headers[1] = "enrg2";
    headers[2] = "enrg4";

    // Stiffness headers
    headers[3] = "stiffx";
    headers[4] = "stiffy";

    return headers;
}



const double* SSEobservables::createValuesArray() {
    static double* values = new double[5];

    // Energy
    values[0] = enrg;
    values[1] = enrg2;
    values[2] = enrg4;

    // Stiffness
    values[3] = stiffx;
    values[4] = stiffy;

    // V_real values (uncomment and fill for the required 12 elements)

    return values;
}




// Function to calculate position vectors with periodic boundary conditions
Site* calculatePositionVectors(int Lx, int Ly) {
    Site* positions = new Site[Lx * Ly];

    // Reference position for site 0
    int refX = 0, refY = 0;

    // Calculate position vectors
    for (int i = 0; i < Lx * Ly; ++i) {
        int relX = i % Lx;
        int relY = i / Lx;
        positions[i].x = relX - refX;
        positions[i].y = relY - refY;

        // Apply periodic boundary conditions along x
        if (positions[i].x > Lx / 2) {
            positions[i].x -= Lx;
        }
        if (positions[i].x < -Lx / 2) {
            positions[i].x += Lx;
        }

        // Apply periodic boundary conditions along y
        if (positions[i].y > Ly / 2) {
            positions[i].y -= Ly;
        }
        if (positions[i].y < -Ly / 2) {
            positions[i].y += Ly;
        }
    }

    return positions;
}

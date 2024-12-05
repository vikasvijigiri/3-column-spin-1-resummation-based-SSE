#include "SSElattice.hpp"
#include "SSEvariables.hpp"
#include "SSEupdates.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>


int** 
SSEvariables::lattice_sites(int lx, int ly, int Nb, int N){
        bsites = new int * [2];
 	for (int i = 0; i < 2; i++) bsites[i] = new int[Ns];
        for (y1=0; y1<ly; ++y1){
            for (x1=0; x1<lz; ++z1){
                s1 = x1+y1*lx;
                x2 = (x1+1) % lx;
                y2 = y1;
                bsites[0][s1] = s1;
                bsites[1][s1] = x2+y2*lx;
                x2 = x1;
                y2 = (y1+1) % ly;
                bsites[0][s1+N] = s1;
                bsites[1][s1+N] = x2+y2*lx;
            }
        }        
        return bsites;
   }     

#include <fstream>
#include <iostream>
#include <iomanip>

#include "../include/SSElattice.hpp"
#include "../include/SSEvariables.hpp"



// constants definition
const int H_BONDS = 4;

// Local Functions declaration
void set_H_Bonds(int** , int , int , int , int , int , int* );
void set_QQ_Bonds(int** , int , int , int , int , int , int , int , int* );


// *************************************************************************
//                              Main 
// *************************************************************************

void
SSEvariables::declare_variables()
{
  // Read variables from file "input_param.dat" 
	std::string INPUT_FILE = "input_params.dat";
	std::ifstream vars(INPUT_FILE);

	vars >> lx;             // Size along x
	vars >> ly;             // Size along y
  vars >> J[0]; 			    // Heisenberg AFM interaction strength
  vars >> J[1];           // QQ AFM interaction strength
  vars >> N;              // N of SU(N).	
  vars >> Beta;           // Inverse temp
	vars >> isteps;         // Equilibration steps
	vars >> iter;           // Measurement steps
 
	
	Ns = lx * ly;				    // No. of spins in a single layer.
	Nm = 2 * Ns;		        // No. of spins totally.
  Np = Ns;

	Nb[0] = 1*4*Ns;			    // No. of Heisenberg bonds
	Nb[1] = 1*8*Ns; 		    // No. of Q_Q bonds

	Lc = 30;                // Initial SSE cutoff length
	n1 = 0;                 // Initial no. of diagonal operators


  // Probabilities to insert and delete
 	prob_in = 0.0;
 	prob_rm = 0.0;
 	cum_prob[0] = (J[0] / N) * (Beta * Nb[0]);    //0.5 * J[0] * Nb[0] *Beta;
 	cum_prob[1] = (J[1] / N) * (Beta * Nb[1]);    //0.5 * 0.5 * J[1] * Nb[1] *Beta;


 	for (int i=0; i<2; i++) prob_in += cum_prob[i];   // insertion probability
 	for (int i=0; i<2; i++) cum_prob[i] /= prob_in;   // removal   probability
 	prob_rm = 1.0/prob_in;	
}


void SSEvariables::print_variables() {
	std::cout << "Lx   : " << lx         << '\n';
	std::cout << "Ly   : " << ly         << '\n';
	std::cout << "JH   : " << J[0]       << '\n';
	std::cout << "JQ   : " << J[1]       << '\n';
	std::cout << "N    : " << N          << '\n';
	std::cout << "Beta   : " << Beta     << '\n';
	std::cout << "Eq steps : " << isteps << '\n';
	std::cout << "Ms steps : " << iter   << '\n';
}

void SSEvariables::lattice_sites()
{

    // stiffness related Info

    JHsgnx = new int[Nb[0]];
    JHsgny = new int[Nb[1]];

    std::fill(JHsgnx, JHsgnx + Nb[0]/2, 1);
    std::fill(JHsgny, JHsgny + Nb[0]/2, 0);
    std::fill(JHsgnx + Nb[0]/2, JHsgnx + Nb[0], 0);
    std::fill(JHsgny + Nb[0]/2, JHsgny + Nb[0], 1);

    JQsgnx = new int[Nb[1]];
    JQsgny = new int[Nb[1]];

    std::fill(JQsgnx, JQsgnx + Nb[1]/2, 1);
    std::fill(JQsgny, JQsgny + Nb[1]/2, 0);
    std::fill(JQsgnx + Nb[1]/2, JQsgnx + Nb[1], 0);
    std::fill(JQsgny + Nb[1]/2, JQsgny + Nb[1], 1);


    // *************************************************

    // Bond Info
    JHsites = new int*[Nb[0]];
    JQsites = new int*[Nb[1]];

    for (int i = 0; i < Nb[0]; ++i) JHsites[i] = new int[2];
    for (int i = 0; i < Nb[1]; ++i) JQsites[i] = new int[4];




    for (int y1 = 0; y1 < ly; ++y1) {
        for (int x1 = 0; x1 < lx; ++x1) {
            // Heisenberg interaction bonds

            // Horizontal bonds
            int s1 = x1 + y1 * lx;
            int x2 = (x1 + 1) % lx;
            int s2 = x2 + y1 * lx;

            int baseIndexHorizontal = 4 * s1;

            int shiftValues13[2] = {0, 0}; // lambda = 1 
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 0, Ns, shiftValues13);//, Hsgn, 1.);            
            
            int shiftValues14[2] = {0, 1}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 1, Ns, shiftValues14);//, Hsgn, lambda);           

            int shiftValues15[2] = {1, 0}; // lambda = 0.5
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 2, Ns, shiftValues15);//, Hsgn, lambda);           

            int shiftValues16[2] = {1, 1}; // lambda = 1
            set_H_Bonds(JHsites, baseIndexHorizontal, s1, s2, 3, Ns, shiftValues16);//, Hsgn, 1.);          


//            // vertical bonds
//            x2 = x1;
//            int y2 = (y1 + 1) % ly;
//            s2 = x2 + y2 * lx;

//            int baseIndexVertical = Nb[0] / 2 + 4 * s1;

//            int shiftValues17[2] = {0, 0}; // lambda = 1
//            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 0, Ns, shiftValues17);//, Hsgn, 1.);              
//            
//            int shiftValues18[2] = {0, 1}; // lambda = 0.5
//            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 1, Ns, shiftValues18);//, Hsgn, lambda);              

//            int shiftValues19[2] = {1, 0}; // lambda = 0.5
//            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 2, Ns, shiftValues19);//, Hsgn, lambda);            

//            int shiftValues20[2] = {1, 1}; // lambda = 1
//            set_H_Bonds(JHsites, baseIndexVertical, s1, s2, 3, Ns, shiftValues20);//, Hsgn, 1.);             


            // three-spin interaction bonds

            // Horizontal bonds
            s1 = x1 + y1 * lx;
            x2 = (x1 + 1) % lx;
            s2 = x2 + y1 * lx;

            int s3 = s2;
            int x3 = (x1 + 2) % lx;
            int s4 = x3 + y1 * lx;

            baseIndexHorizontal = 8 * s1;

            int shiftValues21[4] = {0, 0, 1, 1}; // lambda = 1
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 0, Ns, shiftValues21);//, QQsgn, 1.);           

            int shiftValues22[4] = {1, 1, 0, 0}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 1, Ns, shiftValues22);//, QQsgn, lambda);          

            int shiftValues23[4] = {0, 1, 0, 1}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 2, Ns, shiftValues23);//, QQsgn, lambda);          

            int shiftValues24[4] = {1, 0, 1, 0}; // lambda = 0.25
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 3, Ns, shiftValues24);//, QQsgn, lambda*lambda); 

            int shiftValues25[4] = {0, 1, 0, 0}; // lambda = 1
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 4, Ns, shiftValues25);//, QQsgn, 1.);           

            int shiftValues26[4] = {1, 0, 1, 1}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 5, Ns, shiftValues26);//, QQsgn, lambda);          

            int shiftValues27[4] = {0, 0, 1, 0}; // lambda = 0.5
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 6, Ns, shiftValues27);//, QQsgn, lambda);          

            int shiftValues28[4] = {1, 1, 0, 1}; // lambda = 0.25
            set_QQ_Bonds(JQsites, baseIndexHorizontal, s1, s2, s3, s4, 7, Ns, shiftValues28);//, QQsgn, lambda*lambda);          

//            // Vertical bonds
//            s1 = x1 + y1 * lx;
//            s2 = (s1 + lx) % Ns;

//            s3 = s2;
//            s4 = (s3 + lx) % Ns;

//            baseIndexVertical = Nb[1] / 2 + 8 * s1;

//            int shiftValues31[4] = {0, 0, 1, 1}; // lambda = 1
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 0, Ns, shiftValues31);//, QQsgn, 1.);           

//            int shiftValues32[4] = {1, 1, 0, 0}; // lambda = 0.5
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 1, Ns, shiftValues32);//, QQsgn, lambda);          

//            int shiftValues33[4] = {0, 1, 0, 1}; // lambda = 0.5
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 2, Ns, shiftValues33);//, QQsgn, lambda);          

//            int shiftValues34[4] = {1, 0, 1, 0}; // lambda = 0.25
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 3, Ns, shiftValues34);//, QQsgn, lambda*lambda); 

//            int shiftValues35[4] = {0, 1, 0, 0}; // lambda = 1
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 4, Ns, shiftValues35);//, QQsgn, 1.);           

//            int shiftValues36[4] = {1, 0, 1, 1}; // lambda = 0.5
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 5, Ns, shiftValues36);//, QQsgn, lambda);          

//            int shiftValues37[4] = {0, 0, 1, 0}; // lambda = 0.5
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 6, Ns, shiftValues37);//, QQsgn, lambda);          

//            int shiftValues38[4] = {1, 1, 0, 1}; // lambda = 0.25
//            set_QQ_Bonds(JQsites, baseIndexVertical, s1, s2, s3, s4, 7, Ns, shiftValues38);//, QQsgn, lambda*lambda);           
//   
        }
    }
}


// ****************************************************************************
//                              Local functions 
// ****************************************************************************

void set_H_Bonds(int** bonds, int baseIndex, int s1, int s2, int offset, int N, int* shift) {//, double *sgn, double lmbda)  {
    bonds[baseIndex + offset][0] = s1 + N * shift[0];
    bonds[baseIndex + offset][1] = s2 + N * shift[1];
    //sgn[baseIndex + offset] = lmbda;

}

void set_QQ_Bonds(int** bonds, int baseIndex, int s1, int s2, int s3, int s4, int offset, int N, int* shift) {//, double *sgn, double lmbda) {
    bonds[baseIndex + offset][0] = s1 + N * shift[0];
    bonds[baseIndex + offset][1] = s2 + N * shift[1];
    bonds[baseIndex + offset][2] = s3 + N * shift[2];
    bonds[baseIndex + offset][3] = s4 + N * shift[3];
    //sgn[baseIndex + offset] = lmbda;
}


#include <iostream>
#include <tuple>

#include "../include/SSEvariables.hpp"
#include "../include/SSEobservables.hpp"
#include "../include/SSEupdates.hpp"


// Total effective legs in the system of Heisenberg and 3-spin bonds.
const int H_Tlegs = 8;

static ran rann;
using namespace std;


// ===================================================================
//                            Main
// ===================================================================



// Equilibration
void SSEupdates::Equilibration(SSEvariables& vars) 
{
  auto start = std::chrono::system_clock::now();

	looper(vars);	// Start with a random configuration and wire it.
	for (int i = 0; i < vars.isteps; ++i) EqMCS(vars);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Equilibration in minutes, " << elapsed_seconds.count()/60. << std::endl;
}




// Measurements
void SSEupdates::Measurements(SSEvariables& vars, SSEobservables& obs) 
{
  auto start = std::chrono::system_clock::now();

	obs.Initiate_observables();
	for (int i = 0; i < vars.iter; ++i) MsMCS(vars, obs);

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Measurements in minutes, " << elapsed_seconds.count()/60. << std::endl;
}




// *********************************************************************************************
//                                  Weights part 
// *********************************************************************************************


void
SSEupdates::initvrtx_dirloop(SSEvariables& vars)
{
	vxprb = new int *[4];
	lop = new int [vars.Np];
	for (int i = 0; i < 4; i++)
	{
		vxprb[i] = new int[vars.Np];
	}

	//*************************** Projection operator probabilities *******************************
	for (int ip = 0; ip < vars.Np; ++ip) {		 // diff loop to start with
		vxprb[0][ip] = 3; // Now making it diagonal after prob update.
		vxprb[1][ip] = 2;
		vxprb[2][ip] = 1;
		vxprb[3][ip] = 0;	
		lop[ip] = -1;			// diagonal
	}
}




void
SSEupdates::checkl(SSEvariables& vars)
{
	if (int(vars.n1 * 4 / 3) > vars.Lc)
	{
		int Lc_new = int(vars.n1 * 4 / 3);
		int *str_tmp = new int[Lc_new];

		for (int i = 0; i < vars.Lc; ++i)
		{
				str_tmp[i] = str[i];
		}		
		for (int i = vars.Lc; i < Lc_new; ++i)
		{
				str_tmp[i] = 0;
		}			
		
		delete[] str;
		delete[] X;

		vars.Lc = Lc_new;
		str = new int[vars.Lc];
		X = new int[H_Tlegs*(vars.Lc + vars.Np)];

		for (int i = 0; i < vars.Lc; ++i)
		{
			str[i] = str_tmp[i];
		}
	
		delete[] str_tmp;
	}
	looper(vars);
}


//********************************************************************************************

void
SSEupdates::initialize(SSEvariables& vars)
{
  
  frst  = new int [vars.Nm];
  last  = new int [vars.Nm];
  X     = new int [H_Tlegs*(vars.Lc + vars.Np)];  
  str   = new int [vars.Lc];
  
  
	// H-bond operator part
	for (int j = 0; j < vars.Lc; j++) str[j] = 0;


//  vars.n1 = 5;

//  str[0]=vars.Nb[0] + 0+1;
//  str[2]=vars.Nb[0] + 3+1;
//  str[10]=vars.Nb[0] + 2+1;
//  str[14]=vars.Nb[0] + 1+1;
//  str[17]=vars.Nb[0] + 2+1;

 	type[0] = 1;
	type[1] = 2;
}





// Equilibration Monte Carlo Sweep
void SSEupdates::EqMCS(SSEvariables& vars) 
{
	op_update(vars);
	checkl(vars);
	//update.looper();
} 

// Measurement Monte Carlo Sweep
void SSEupdates::MsMCS(SSEvariables& vars, SSEobservables& obs) 
{
	op_update(vars);
	obs.measure_observables(vars);
} 



//***********************************************************************************************



/*
int 
SSEupdates::insert_op_dnl (int v0, int v1)
{
  //std::cout << v0 << "  " << v1 << std::endl;
  
  // ******************************************* Check SAME loop or NOT. ***********************************************
  int vi = v0;
  int l,il;
  int ic, oc;
  double r, p;
  int dnl = 0;
  
 	while (1) {	  		
      ic = vi % 4;
  		l = int(vi / 4);   
  		//std::cout << vi << std::endl;    
      if (l < Lc) {  // When loop encounters Heisenberg bonds.
        oc = ic ^ 1;
      } else {
      	oc = vxprb[ic][l-Lc];
      }
      vi = 4 * (int(vi / 4)) + oc;
  		if (vi == v1)	// same loop 
  		{
  			dnl = 1;
  			break;
  		}		
  		vi = X[vi];
  		//std::cout << vi << std::endl;    
  		if (vi == v0)	// diff loop
  		{
  			dnl = -1;
  			break;
  		}
  		if (vi == v1)	// same loop 
  		{
  			dnl = 1;
  			break;
  		}		
  }
  //std::cout << " dnl is "<< dnl << std::endl;
	return dnl;
}
*/




	/*
	int b, o, s1, s2, dnl;
	double p ;
	for (int i = 0; i < Lc; ++i)
	{
		o = str[i];
		if (o == 0)
		{
			b = int(rann() *(Nb[0]));
			s1 = JHsites[b][0];
			s2 = JHsites[b][1];
			auto [v1, v2, pos_type1, pos_type2] = get_legs(s1, s2, i, true);
			update_links(s1, s2, i, v1, v2, pos_type1, pos_type2, true);
 			dnl = calculate_dnl(i,3);
     	p = prob_in * pow(N, -dnl) / float(Lc - n1);
      if (p > rann()){
		     		str[i] = b + 1;
		     		str2[i] = 0;
		      	n1 += 1; 
		  }
		  else {
		  			update_links(s1, s2, i, v1, v2, pos_type1, pos_type2, false);
		  }
		}
		else 
		{
			b  = str[i] - 1;
	    s1 = JHsites[b][0];
	    s2 = JHsites[b][1];
			auto [v1, v2, pos_type1, pos_type2] = get_legs(s1, s2, i, false);
 			dnl = calculate_dnl(i,3);			
			p = prob_rm * pow(N, dnl) * float(Lc-n1+1);
			if (p > rann()){
					str[i] = 0;
					str2[i] = 0;
					n1 -=  1; 
					update_links(s1, s2, i, v1, v2, pos_type1, pos_type2, false);
			}
    }
  }

  // Projector update (or) rewiring.
	for (int i = 0; i < Np; ++i)
	{
 			int dnl = calculate_dnl(i+Lc,1);
     	p =  pow(N, dnl);
      if (p > rann()) {rewire(i);}
	}	
	*/
	



    /*
    if (o > 0) {
      b = str[i] - 1;
      s1 = bsites[b][0];
      s2 = bsites[b][1];

      v1 = last[s1];
      v2 = last[s2];
      //lpos[int(v0 / 4)] = i;
      if (v1 != -1) {
        X[v1] = v0;
        X[v0] = v1;
      } else {
        frst[s1] = v0;
      }
      if (v2 != -1) {
        X[v2] = v0 + 1;
        X[v0 + 1] = v2;
      } else {
        frst[s2] = v0 + 1;
      }
      last[s1] = v0 + 2;
      last[s2] = v0 + 3;
    }
    */

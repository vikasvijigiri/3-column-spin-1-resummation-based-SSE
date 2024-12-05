#include <iostream>
#include <tuple>

#include "../include/SSEvariables.hpp"
#include "../include/SSEupdates.hpp"




// Total effective legs in the system of Heisenberg and 3-spin bonds.
const int H_Tlegs = 8;


// =======================================================
//                      Main 
// =======================================================

void
SSEupdates::looper(SSEvariables& vars)
{
	// Linked-list storage 
	initialize_links(vars);
  link_interaction_bonds(vars);
  link_projector_bonds(vars);
  link_pbc_bonds(vars);
  //print_links(vars);

}



// =======================================================
//                   Local Functions 
// =======================================================

void SSEupdates::initialize_links(SSEvariables& vars) {
	std::fill_n(frst, vars.Nm, -1);
	std::fill_n(last, vars.Nm, -1);
	std::fill_n(X, H_Tlegs * (vars.Lc + vars.Np), -1);
}


void SSEupdates::print_links(SSEvariables& vars) {
	for (int i=0; i < H_Tlegs*(vars.Lc + vars.Np); i++)
		std::cout << i << "  <--> " << X[i] << std::endl;
}


void SSEupdates::link_pbc_bonds(SSEvariables& vars) {
	for (int k = 0; k < vars.Nm; ++k) {
		int v1 = frst[k];
		if (v1 != -1) {
			int v2 = last[k];
			X[v2] = v1;
			X[v1] = v2;
		}
	}
}


void SSEupdates::link_interaction_bonds(SSEvariables& vars) {
  int  s1, s2;
  for (int i = 0; i < vars.Lc; ++i) {
    int v0 = H_Tlegs*i;
    int o  = str[i];
    if (o != 0) {
      int b = (o - 1 < vars.Nb[0]) ? o - 1 : (o - 1 - vars.Nb[0]); //(o - 1) % vars.Nb[0];
      int v = (o - 1 < vars.Nb[0]) ? 0 : 1;   //(o - 1) / vars.Nb[0];//str2[i];
      for (int j = 0; j < type[v]; j++) {
        if (v == 0) {
          s1 = vars.JHsites[b][2 * j];
          s2 = vars.JHsites[b][2 * j + 1];
        } else if (v == 1) {
          s1 = vars.JQsites[b][2 * j];
          s2 = vars.JQsites[b][2 * j + 1];
        } else {
          printf("wrong vertex ");
          exit(0);
        }
        //std::cout << s1 << "   " << s2 << std::endl;
        int v1 = last[s1];
        int v2 = last[s2];
        if (v1 != -1) {
          X[v1] = v0 + 4*j;
          X[v0 + 4*j] = v1;
        } else {
          frst[s1] = v0 + 4*j;
        }
        if (v2 != -1) {
          X[v2] = v0 + 1 + 4*j;
          X[v0 + 1 + 4*j] = v2;
        } else {
          frst[s2] = v0 + 1 + 4*j;
        }
        last[s1] = v0 + 2 + 4*j;
        last[s2] = v0 + 3 + 4*j;
      }
    }    
  }	
}


void SSEupdates::link_projector_bonds(SSEvariables& vars) {
	for (int i = 0; i < vars.Np; ++i) {  //Projector-bonds
	  int v0 = H_Tlegs*(vars.Lc + i);
	  int s1 = i;
	  int s2 = i + vars.Np;
	  int v1 = last[s1];
	  int v2 = last[s2];

	  if (v1 != -1) {
		  X[v1] = v0;
		  X[v0] = v1;
	  } else {
		  frst[s1] = v0;
	  }
	  if (v2 != -1) {
		  X[v2] = v0 + 1;
		  X[v0 + 1] = v2;
	  }
	  else{
		  frst[s2] = v0 + 1;
	  }
	  last[s1] = v0 + 2;
	  last[s2] = v0 + 3;
	}
}


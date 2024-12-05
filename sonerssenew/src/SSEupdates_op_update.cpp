#include <iostream>
#include <tuple>

#include "../include/SSEvariables.hpp"
#include "../include/SSEupdates.hpp"


// Total effective legs in the system of Heisenberg and 3-spin bonds.
const int H_Tlegs = 8;

static ran rann;
using namespace std;






// ======================================================
//                   Main 
// ======================================================

void
SSEupdates::op_update(SSEvariables& vars)
{
	// Bond operators update. 
  for (int i = 0; i < vars.Lc; ++i) {
    int o = str[i]; 
    //std::cout << i << "  " << vars.n1 << std::endl;
    if (o == 0) { // Encountered Identity, try to insert the interaction bond operators
        //std::cout << "A" << std::endl;        
        insert_operator(vars, i);
        //std::cout << "B" << std::endl;
    } 
    else  {
        remove_operator(vars, i);
    } 
  }

  // Projector update (or) rewiring.
	for (int i = 0; i < vars.Np; ++i){
		int dnl = calculate_dnl(vars.Lc, i + vars.Lc, 0, 1);
   	double p =  pow(vars.N, dnl);
    if (p > rann()) rewire(i);
	}	
}



// ======================================================
//                   Local functions 
// ======================================================


void SSEupdates::insert_operator(SSEvariables& vars, int i) {
  double r = rann();
  double cp = 0.0;
	int ss1[2], ss2[2], vl1[2], vl2[2], pos_typel1[2], pos_typel2[2];
  bool false_entered[2] = {false, false};

  // Insert the operator
  for (int k = 0; k < 2; k++) {
    cp += vars.cum_prob[k];
    if (cp > r and vars.cum_prob[k] > 1e-6) { // J or Q or Biquad?
      int b = int(rann() * vars.Nb[k]);
      for (int j = 0; j < type[k]; ++j) {
        double p = 1.0;
        if (k == 0) {
          ss1[j] = vars.JHsites[b][2 * j];
          ss2[j] = vars.JHsites[b][2 * j + 1];
        } else if (k == 1) {
          ss1[j] = vars.JQsites[b][2 * j];
          ss2[j] = vars.JQsites[b][2 * j + 1];
        }
        ////std::cout << i<< "    "  << ss1[j] << "  " << ss2[j] << "  " << i << "  " << j << std::endl;
			  auto [v1, v2, pos_type1, pos_type2] = get_legs(vars.Lc, vars.Np, ss1[j], ss2[j], i, j, true);
			  vl1[j] = v1;
			  vl2[j] = v2;
			  pos_typel2[j] = pos_type2;
			  pos_typel1[j] = pos_type1;
			  update_links(ss1[j], ss2[j], i, j, v1, v2, pos_type1, pos_type2, true); // fake update
        int dnl = calculate_dnl(vars.Lc, i, j, 3);
 				p = p * pow(vars.N, -dnl);
        p = p * vars.prob_in / double(vars.Lc - vars.n1);
        if (p < rann()) { 
          update_links(ss1[j], ss2[j], i, j, v1, v2, pos_type1, pos_type2, false); 
          false_entered[j] = true;
        }
      }
      if (false_entered[1] and !false_entered[0]) {
          update_links(ss1[0], ss2[0], i, 0, vl1[0], vl2[0], pos_typel1[0], pos_typel2[0], false);
      }
//      if (!false_entered[1] and false_entered[0]) {
//          update_links(ss1[1], ss2[1], i, 1, vl1[1], vl2[1], pos_typel1[1], pos_typel2[1], false);
//      }
      if (!false_entered[1] and !false_entered[0]) {vars.n1 += 1; str[i] = vars.Nb[0]*k + b+ 1;  }
    break;
    }
  }
}


void SSEupdates::remove_operator(SSEvariables& vars, int i) {
  int o = str[i] - 1;  
  int b = (o < vars.Nb[0]) ? o : (o - vars.Nb[0]);
  int k = (o < vars.Nb[0]) ? 0 : 1;

	int ss1[2], ss2[2], vl1[2], vl2[2], pos_typel1[2], pos_typel2[2];
  bool false_entered[2] = {false, false};

  // Remove the operator
  for (int j = 0; j < type[k]; ++j) {
    double p = 1.0;
    if (k == 0) {
      ss1[j] = vars.JHsites[b][2 * j];
      ss2[j] = vars.JHsites[b][2 * j + 1];
    } else if (k == 1) {
      ss1[j] = vars.JQsites[b][2 * j];
      ss2[j] = vars.JQsites[b][2 * j + 1];
    }
    ////std::cout << ss1[j] << "  " << ss2[j] << "  " << i << "  " << j << std::endl;
	  auto [v1, v2, pos_type1, pos_type2] = get_legs(vars.Lc, vars.Np, ss1[j], ss2[j], i, j, false);
	  vl1[j] = v1;
	  vl2[j] = v2;
	  pos_typel2[j] = pos_type2;
	  pos_typel1[j] = pos_type1;						
	  int dnl = calculate_dnl(vars.Lc, i, j, 3);
		p = p * pow(vars.N, dnl);
    p = p * vars.prob_rm * double(vars.Lc - vars.n1 + 1);
    if (rann() < p) {update_links(ss1[j], ss2[j], i, j, v1, v2, pos_type1, pos_type2, false); false_entered[j] = true; 
    if (k==0) false_entered[1]=true;}
  }
  if (false_entered[1] and !false_entered[0]) update_links(ss1[0], ss2[0], i, 0, vl1[0], vl2[0], pos_typel1[0], pos_typel2[0], true);
  if (!false_entered[1] and false_entered[0]) update_links(ss1[1], ss2[1], i, 1, vl1[1], vl2[1], pos_typel1[1], pos_typel2[1], true);
  if (false_entered[0] and false_entered[1]) {vars.n1 -= 1; str[i] = 0;}
}



std::tuple <int, int, int, int> SSEupdates::get_legs (int Lcd, int Npd, int s1, int s2, int i, int i1, bool insert) {

  
  int v[2], v1[2], w[2], pos_type[2], l0, l2, id=0, iu=0;

  v[0] = frst[s1];  // Below
  v[1] = frst[s2];
  
  w[0] = last[s1];	// Above
  w[1] = last[s2];
  
  v1[0] = v[0];
  v1[1] = v[1];
  
  l0 = H_Tlegs*i + 4*i1;
  l2 = l0 + 2;

  
  // ********************************** Getting the leg vertex of the site. **************************************
  if (insert) {
    for (int j = 0; j < 2; j++) {
  		if ( i < int(v[j]/H_Tlegs) and i < int(w[j]/H_Tlegs) )	// link ABOVE to the current one
  			{	
  				v1[j] = v[j];
  				pos_type[j] = 0;
  			}
  			else if ( i > int(v[j]/H_Tlegs) and i < int(w[j]/H_Tlegs) )	// link SANDWITCH to the current one.
  			{
  				v1[j] = v[j];
  				while (1)
  				{
  					v1[j] += 2;
  					v1[j] %= H_Tlegs*(Lcd + Npd); 
  					v1[j] = X[v1[j]];
  					if (int(v1[j]/H_Tlegs) > i)
  					{
  						v1[j] = X[v1[j]];
  						break;
  					}
  				}
  				pos_type[j] = 1;
  			}
  			else { // Else terminate with error msg.
    			cout << "INSERT: legs are   " << j << "  " << i << "  " << int(v[j]/H_Tlegs) << "  " << int(w[j]/H_Tlegs) << endl;   		 			
  	      printf(" Error in getting vertex leg! ");
					exit(0);  			
  			}
   	}
  }
  else {
   for (int j = 0; j < 2; j++)
  	{
         	id = int(X[(l0+j)]/H_Tlegs);
        	iu = int(X[(l2+j)]/H_Tlegs);  	
					if ( (i > id and i < iu) or (i < id and i > iu) ){
						pos_type[j] = 1; // Sandwich
					}
					else if ( i < id and i < iu ){
						pos_type[j] = 0; // Above		
					}
					else { // Else terminate with error msg.
						cout << "REMOVAL: legs are   " << i << "  " << v1[0] << "  " << v1[1] << "  " << id << "  " << iu << endl;   		
						printf(" Error in computing vertex leg! ");
						exit(0);  			
					}				
					v1[j] = X[l0 + j];
		}  		
  }
  return {v1[0], v1[1], pos_type[0], pos_type[1]};		
}





int SSEupdates::calculate_dnl (int Lcd, int i, int i1, int shift)
{
	int v0 = H_Tlegs*i + 4*i1;
	int v1 = v0 + shift;
	
  int vi = v0;
  int oc;
  int dnl = 0;
  int p = 0;

  while (true) {		
      int ic = vi % 4;
      int l = int(vi / H_Tlegs);
      if (l < Lcd) {  // When loop encounters Heisenberg bonds.
        oc = ic ^ 1;
        p = 1;
      } else {
      	oc = vxprb[ic][l-Lcd];
      	p = 0;
      }
     
      vi = H_Tlegs * (int(vi / H_Tlegs)) + oc +  p*4*int((vi % H_Tlegs)/4);  
  		vi = X[vi];        
  		if (vi == v1)	// same loop 
  		{
  			dnl = 1;
  			break;
  		}	  		      	             
  		if (vi == v0)	// diff loop
  		{
  			dnl = -1;
  			break;
  		}  		
  }
  return dnl;	
}


void SSEupdates::update_links (int s1, int s2, int i, int i1, int vl1, int vl2, int pos_type1, int pos_type2, bool insert)
{
  int s[2], v1[2], vu[2], vd[2], pos_type[2];
  int l0, l2;
  
  s[0] = s1;
  s[1] = s2;
  
  v1[0] = vl1;
  v1[1] = vl2;

  pos_type[0]=pos_type1;
  pos_type[1]=pos_type2; 
  
    
  l0=H_Tlegs*i+4*i1;
  l2=l0+2;

  for (int j = 0; j < 2; j++){
	  vd[j] = v1[j];
	  if (pos_type[j] == 0)	    // link above.  
	  {
	    if (insert){
	      vu[j] = X[vd[j]];	    

	      X[l0+j]=vu[j];
	      X[vu[j]]=l0+j;
	        
	      X[l2+j]=vd[j];
	      X[vd[j]]=l2+j;
	    
	      frst[s[j]]=l0+j;
	    	//cout << "links are above ins  "  << vd[j] << "  " << vu[j] << "  " << X[vd[j]] << "  " << X[vu[j]] << "  " << X[X[vd[j]]] << "  " << X[X[vu[j]]] 
	    	// << endl;     	    	        	        
	    }
	    else {
				vd[0] = X[l0 + j];
				vd[1] = X[l2 + j];

				frst[s[j]]=vd[1];

				X[vd[0]] = vd[1];
				X[vd[1]] = vd[0];
							
				X[l0 + j] = -1;
				X[l2 + j] = -1;
	    	//cout << "links are above rem  "  << vd[0] << "  " << vd[1] << "  " << X[vd[0]] << "  " << X[vd[1]] << "  " << X[X[vd[0]]] << "  " << X[X[vd[1]]] 
	    	//<< endl;     	    	        	        
	    }
	    //cout << "links are above    " << X[29] << "  " << X[31] << "  " << frst[s1] << endl;     	    	        	    	    
	  }
	  else if (pos_type[j] == 1) 	    // link Sandwitching.
	  {
	    if (insert){	  
	      vu[j] = X[vd[j]];
	    
	      X[l0+j]=vd[j];
	      X[vd[j]]=l0+j;
	    
	      X[l2+j]=vu[j];
	      X[vu[j]]=l2+j;
	    	//cout << "links are sandw ins  "  << vd[j] << "  " << vu[j] << "  " << X[vd[j]] << "  " << X[vu[j]] << "  " << X[X[vd[j]]] << "  " << X[X[vu[j]]] 
	    	// << endl;     	    	        
	    }
	    else {
				vd[0] = X[l0 + j];
				vd[1] = X[l2 + j];

				X[vd[0]] = vd[1];
				X[vd[1]] = vd[0];
					
				X[l0 + j] = -1;
				X[l2 + j] = -1;
	    	//cout << "links are sandw rem  "  << vd[0] << "  " << vd[1] << "  " << X[vd[0]] << "  " << X[vd[1]] << "  " << X[X[vd[0]]] << "  " << X[X[vd[1]]] 
	    	//<< endl;      	    	        	        
	    } 
	  }  
	  else {
	      printf(" Error in the pos_type index! ");
				exit(0);	  
	  }	  
  }  
}


void SSEupdates::rewire(int ip)
{
	// Update...
	//*************************** Projection operator probabilities *******************************
	if (lop[ip] == -1){ // it was diagonal
		vxprb[0][ip] = 2; // Now making straight after prob update.
		vxprb[1][ip] = 3;
		vxprb[2][ip] = 0;
		vxprb[3][ip] = 1;		
		lop[ip] = 1;			// straight
	} else { 				  		// it was straight
		vxprb[0][ip] = 3; // Now making it diagonal after prob update.
		vxprb[1][ip] = 2;
		vxprb[2][ip] = 1;
		vxprb[3][ip] = 0;	
		lop[ip] = -1;			// diagonal
	}
}





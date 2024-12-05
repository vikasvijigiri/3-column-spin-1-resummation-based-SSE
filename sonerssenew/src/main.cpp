#include <iostream>
//#include <chrono>

//#include "SSElattice.hpp"
#include "../include/SSEvariables.hpp"
#include "../include/SSEobservables.hpp"
#include "../include/SSEupdates.hpp"			//Headers
#include "../include/writeresults.hpp"


//static ran rann;


int main()
{
  // Creating all necessary instances
	SSEvariables vars;	
	SSEupdates update;
	SSEobservables obs;
	writeresults write;

  // Declare variables;
	vars.declare_variables();
  vars.print_variables();
  vars.lattice_sites();

  // Initialize directed loop weights for mini-spins space (Projection part)
	update.initvrtx_dirloop(vars);
	update.initialize(vars);
	//update.weights();		 	

  // Main QMC algorithms
  update.Equilibration(vars);
  update.Measurements(vars, obs);

  // write output to a file
  write.save_bin_data("RSSE_data.txt", obs.createHeadersArray(), obs.createValuesArray(), 5);

	return 0;
}

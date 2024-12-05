#ifndef _SSEUPDATES_HPP_DEFINED_
#define _SSEUPDATES_HPP_DEFINED_

#include "../include/SSEobservables.hpp"

class SSEupdates {
	public: 
	SSEupdates(){}

	int *str;     //, *str2;//, *pstr;
	int type[2];
 	int **vxprb;
	int *lop; 
	int *frst, *last, *X;
	
	void initvrtx_dirloop(SSEvariables& );
	void checkl(SSEvariables& );
	void op_update(SSEvariables& );
	void looper(SSEvariables& );		
  void initialize(SSEvariables& );

  // define Local functions
  void initialize_links(SSEvariables& );
  void link_interaction_bonds(SSEvariables& );
  void link_projector_bonds(SSEvariables& );
  void link_pbc_bonds(SSEvariables& );
  void print_links(SSEvariables& );

  void insert_operator(SSEvariables& , int);
  void remove_operator(SSEvariables& , int);

  std::tuple <int, int, int, int> get_legs (int, int, int, int, int, int, bool );	 
  int  calculate_dnl (int, int, int, int);	
  void update_links (int, int, int, int, int, int, int, int, bool );
  void rewire (int ); 

  void EqMCS(SSEvariables& vars);
  void MsMCS(SSEvariables& vars, SSEobservables& );

  // Global funcs
  void Equilibration(SSEvariables& );
  void Measurements(SSEvariables& , SSEobservables& );
      	
};
#endif


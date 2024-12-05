#ifndef _SSELATTICE_HPP_DEFINED_
#define _SSELATTICE_HPP_DEFINED_
#pragma once
class SSElattice {
  	private:
    	int spin;
  	public:
    	SSElattice() {spin = 1;}
  	SSElattice(int s) {spin = s;}
  	int S(void) {return spin;}
  	void set_S(int s) {spin = s;}
  	void flip(void) {spin *= -1;}
};
#endif

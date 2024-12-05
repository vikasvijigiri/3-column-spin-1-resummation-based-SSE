#ifndef _SSEVARIABLES_HPP_DEFINED_
#define _SSEVARIABLES_HPP_DEFINED_

#include <random>
#include <chrono>

using namespace std::chrono;

#pragma once
class SSEvariables {
  public:
    int lx, ly, Ns, Nm, Np, Nb[2];
    double J[2], N;
    int n1, Lc, isteps, iter;
    double Beta;
    int **JHsites, **JQsites;
    int *JHsgnx, *JHsgny, *JQsgnx, *JQsgny;
	  double prob_in, prob_rm, cum_prob[2]; 	


  	void lattice_sites();
  	void declare_variables();
    void print_variables();


};

class ran {
  private:
    std::mt19937 mt;
    std::uniform_real_distribution < double > dist;

  public:
    ran(double lower = 0.0, double upper = 1.0): mt(std::chrono::high_resolution_clock::now().time_since_epoch().count()), dist(lower, upper){}
    double operator()() {
    return dist(mt);
  }
};
#endif

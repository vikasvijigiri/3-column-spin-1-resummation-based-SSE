#ifndef _SSEOBSERVABLES_HPP_DEFINED_
#define _SSEOBSERVABLES_HPP_DEFINED_


#pragma once
class SSEobservables {
  	public:

			  // Observables

        // Energy
        double enrg;
        double enrg2;
        double enrg4;

        // Stiffness
        double stiffx;
        double stiffy;

        // Observable functions
        // void observables_processing(SSEvariables& );
        void Initiate_observables();
        void measure_observables(SSEvariables& ); 
        const char** createHeadersArray();
        const double* createValuesArray();


};
#endif

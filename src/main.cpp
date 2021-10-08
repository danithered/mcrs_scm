#include <iostream>
//#include <ctime>
#include "randomgen.h"
#include "ca.h"
#include "parameters.h"
#include "rnarep.h"

using namespace std;


/* return values:
 0: Ok
 1: died
 -1: error in reading in argoments
 -2: couldnt open output 
 -3: rng initialisation
*/
int main(int argc, char *argv[]) {
	//Argoments
	if ( Args(argc, argv) ) {
		return(-1);
	}

	//initialise rng
	time_t timer;
	if(std::strlen(par_seed_file)){
		if(randomszam_olvasas(par_seed_file, r)) return -3;
	}
	else { // no seed file specified
		if(par_seed < 0){ // init with time
			randomszam_inic(time(&timer) + par_seed_plus, r);
		}
		else randomszam_inic(par_seed, r); // init with exact seed
	}

	//report init
	cout << "Starting to init simulation " << par_ID << " at " << ctime(&timer); 

	//start to do stuff
	cadv::CellAut automata(par_nrow, par_ncol); //initialise automata

	rnarep::CellContent::patterns.readFile(par_str_pool); //read in pattern file

	automata.neighInic(MARGOLUS_NEIGH, cadv::torus, 0); //init diffusional neighbourhood for Toffoli-Margoulus algorithm
	automata.neighInic(par_Nmet, cadv::torus, 1); //init metabolic neighbourhood 
	automata.neighInic(par_Nrep, cadv::torus, 2); //init replication neighbourhood

	//load if needed
	if(std::strlen(par_load) > 0) automata.init_fromfile(par_load);

	//open output
	if(automata.openOutputs()) { //returns not 0 if fails
		gsl_rng_free(r);

		//report closing
		timer = time(0);
		std::cout << "Simulation " << par_ID << " ending at: " << ctime(&timer) << std::endl << "It had init problems." << std::endl;

		return -2;
	}
	
	//save parameters
	std::string paramfilename(automata.savedir.c_str());
	paramfilename += "/parameters.txt";
	paramsToFile(paramfilename.c_str());

	//Running simulation
	if (automata.rUpdate(par_maxtime)){
		//close rng
		gsl_rng_free(r);

		//report closing
		timer = time(0);
		std::cout << "Simulation " << par_ID << " ending at: " << ctime(&timer) << std::endl << "It has survived." << std::endl;

		return 0; // it has survived
	}
	else{	
		// died out
		//close rng
		gsl_rng_free(r);

		//report closing
		timer = time(0);
		std::cout << "Simulation " << par_ID << " ending at: " << ctime(&timer) << std::endl << "It has died out." << std::endl;

		return 1;
	}

} 

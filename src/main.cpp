#include <iostream>
//#include <ctime>
#include "randomgen.h"
#include "cm.h"
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
	cmdv::CompartPool automata(par_poolsize); //initialise automata

	rnarep::CellContent::patterns.readFile(par_str_pool); //read in pattern file

	//automata.neighInic(MARGOLUS_NEIGH, cadv::torus, 0); //init diffusional neighbourhood for Toffoli-Margoulus algorithm
	//automata.neighInic(par_Nmet, cadv::torus, 1); //init metabolic neighbourhood 
	//automata.neighInic(par_Nrep, cadv::torus, 2); //init replication neighbourhood

	//load if needed
	if(std::strlen(par_load) > 0) automata.init_fromfile(par_load);
	if(std::strlen(par_bubbles) > 0) automata.discoverComparts(par_bubbles);

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
	int endstate = automata.oUpdate(par_maxtime);

	//close rng
	gsl_rng_free(r);

	//report closing
	timer = time(0);
	std::cout << "Simulation " << par_ID << " ending at: " << ctime(&timer) << std::endl; 

	switch(endstate){
		case 0:
			std::cout << "It has survived." << std::endl;
			return 0; // it has survived
		case 1:
			std::cout << "It has died out." << std::endl;
			return 1;
		default:
			if(rnarep::CellContent::no_replicators) std::cout << "It has survived.";
			else std::cout << "It has died out.";
			std::cout << " And quitted with condition: " << endstate << std::endl;
			return endstate;
	}

	return -100;

}


#include <iostream>
#include "randomgen.h"
#include "cm.h"
#include "parameters.h"
#include "rnarep.h"

using namespace std;

int main(int argc, char *argv[]) {

	// INIT
	par_substitution = par_insertion = par_deletion = 0;
	strcpy(par_str_pool, argv[1]);

	randomszam_inic(1, r); // init with exact seed
				      
	rnarep::CellContent::patterns.readFile(argv[1]); //read in pattern file
	cmdv::CompartPool automata(1); //initialise automata

	// LOOK TRU FILES
	for(unsigned int file = 2; file < argc; ++file){
		if(!automata.compartFromFile(argv[file])) continue; // import compart
		
		cmdv::Compart  *comp = *(automata.comparts);

		if( comp->alive() ){
			double sum = 0.0, sumsquared = 0.0;

			for(auto &rep : comp->reps){
				sum += rep->getR();
				sumsquared += std::pow(rep->getR(),2);
			}

			std::cout << argv[file] << "\tM=" << comp->get_M() << "\tSD(R)=" << dvtools::sd(comp->reps.size(), sum, sumsquared) << std::endl;
		} 
//		else std::cout << argv[file] << " is dead" << std::endl;
	}

	// END
	
	//close rng
	gsl_rng_free(r);

	return 0;
}


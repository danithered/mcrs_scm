#include <iostream>
#include "parameters.h"
#include "rnarep.h"
#include "cm.h"
#include "randomgen.h"

using namespace std;

int main(int argc, char *argv[]) {
	//Argoments
	//if ( Args(argc, argv) ) {
	//	return(-1);
	//}

	par_substitution = par_insertion = par_deletion = 0;
	strcpy(par_str_pool, argv[1]);

	randomszam_inic(1, r); // init with exact seed
				      
	rnarep::CellContent::patterns.readFile(argv[1]); //read in pattern file
	cmdv::CompartPool automata(1); //initialise automata



	for(unsigned int file = 2; file < argc; ++file){
		if(!automata.compartFromFile(argv[file])) continue; // import compart
		
		if( (**(automata.comparts)).alive() ){
			std::cout << argv[file] << std::endl;
		} 
		else std::cout << argv[file] << " is dead" << std::endl;
	}

	//close rng
	gsl_rng_free(r);

	return 0;
}


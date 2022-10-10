#include <iostream>
#include "ca.h"
#include "rnarep.h"
#include "randomgen.h"

using namespace std;

/*argoments:
nrow
ncol
noEA
pattern pool
grid file
*/

//int par_insertion, par_deletion, par_substitution;

int main(int argc, char *argv[]) {
	randomszam_inic(200, r); // have to put it here because CellAutt:init calls for it

	par_insertion = par_deletion = par_substitution = 0;
	int par_nrow = std::atoi(argv[1]), par_ncol = std::atoi(argv[2]), par_noEA = std::atoi(argv[3]);
	std::string emptystring("N\t0\t0\t0\t0\t0\t0\t-1");

	//create remp cellcontent file
	rnarep::CellContent cell;

	//load grid int automata
	cadv::CellAut automata(par_nrow, par_ncol); //initialise automata
	rnarep::CellContent::patterns.readFile( argv[4] ); //read in pattern file
	automata.init_fromfile( argv[5] ); //load from file

	//finish emptystring with trailing 0-s for missing enz acts
	for(int aa = par_noEA; aa--;) emptystring += "\t0";

	//running thru grid
	for(int i=0; i < par_nrow*par_ncol; i++){
		cell.die();
		cell.replicate(  *(automata.get(i)->vals)   );
		
		//std::cout << cell.get_seq() << cell.get_str() << 
		if(cell.empty) cout << emptystring << std::endl; 
		else {
			std::cout << *(cell.get_seq())
				<< '\t' << cell.get_str()
				<< '\t' << cell.get_mfe() 
				<< '\t' << cell.getPfold()
				<< '\t' << cell.Pdeg 
				<< '\t' << cell.get_no_sites()
				<< '\t' << cell.getR()
				//<< '\t' << cell.M()
				<< '\t' << cell.get_type() ; 
			for(double *a = cell.geta(), *a_until = cell.geta() + par_noEA; a != a_until; a++){
				cout << '\t' << *a;
			}

			cout << std::endl;
		}
	}


	gsl_rng_free(r);

	return 0;
}

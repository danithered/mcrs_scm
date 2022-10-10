#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
}
#include "randomgen.h"

#include <vector>
#include "rnarep.h"
#include "annot.h"
#include "parameters.h"
#include "dv_tools.h"

using namespace std;

char bases[] = "AGCU"; 
/*
int main(int argc, char *argv[]){
	int lambda = 45, no_seqs=5000000;
	char seq_file[] = "OUT/randseqs.txt";


	//initialise rng
	time_t timer;
	randomszam_inic( time(&timer) , r);

	//open seq file
	std::ofstream file(seq_file, std::ios_base::app);
	if(!file.is_open()) std::cerr << "ERROR: file (" << seq_file << ") not found!" << std::endl;
	
	//intit vars
	std::string seq;
	double mfe;
	char *str;
	str = (char *) vrna_alloc(sizeof(char) *  200 );

	//start computing
	for(; no_seqs--; ){
		//get a sequence
		seq.clear();
		for(int rb=gsl_ran_poisson(r, lambda); rb--;) seq.push_back(bases[gsl_rng_uniform_int(r, 4)]);
		
		//compute str and mfe
		mfe = vrna_fold(seq.c_str(), str);
		file  << seq << '\t' << str << '\t' << seq.length() << '\t' << mfe << std::endl; 
	}






	file.close();
	delete [] (str);

	//close rng
	gsl_rng_free(r);
	return 0;
}
*/


int main(int argc, char *argv[]){
	//Argoments
	if ( Args(argc, argv) ) {
		return(-1);
	}

	par_insertion = par_deletion = par_substitution = 0;

	//initialise rng
	time_t timer;
	randomszam_inic( time(&timer) , r);






	//stuff
	std::string seq_file_start("IN/str/"), seq_file, compl_seq_file, seq;
	rnarep::CellContent::patterns.readFile(par_str_pool); //read in pattern file
	rnarep::CellContent replicator, compl_replicator;
	int no_seqs = 100000000, lambda = 45;
	//int no_seqs = 100000000, lambda = 30;
	//int no_seqs = 100, lambda = 45;
	unsigned long long int type;

//	rnarep::CellContent::patterns.printRules();

	//create output file streams
	std::vector<std::ofstream> files, compl_files;
	std::ofstream tempof;
	double no_types = std::pow(2.0, (double) par_noEA);
	for(int of = 0; of < par_noEA; of++){ 
		seq_file.erase();
		seq_file = seq_file_start;
		seq_file += std::to_string(of+1) ;
		seq_file += "/randseqs_ea" ;
		seq_file += std::to_string(of);
		compl_seq_file = seq_file;
		seq_file += ".txt";
		compl_seq_file += "_compl.txt";

		tempof.open(seq_file, std::ios_base::app);
		files.push_back(std::move(tempof));

		tempof.open(compl_seq_file, std::ios_base::app);
		compl_files.push_back(std::move(tempof));
//		files[of] << "vmi";
		//if(!files[of].is_open()) std::cerr << "ERROR: file (" << seq_file << ") not found!" << std::endl;
	}

	for(; no_seqs--; ){
		//get a sequence
		seq.clear();
//		std::cout << "seq cleared" << seq << std::endl;
		for(int rb=gsl_ran_poisson(r, lambda); rb--;) {
//			std::cout << "add base" << std::endl;
			seq.push_back(bases[gsl_rng_uniform_int(r, 4)]);
		}
//		std::cout << "seq is " << seq << std::endl;
		replicator = seq;

//		std::cout << replicator.getR() <<  std::endl;

		//put it in file
		type = replicator.get_type();
//		std::cout << type << " " << std::log2(type) << seq << std::endl;
//		if(type < no_types) { //dont forget to delete this condition! it just slows it down...
			//if(type && (type % 2 == 0) ) files[(int) std::log2(type) ] << seq << '\t' << replicator.get_str() << '\t' << seq.length() << '\t' << replicator.get_mfe() << std::endl; //this outputs PURE enzimes (so not promiscous ones)

			//this outputs any replicator, which has given activity (note that if it is a promiscous one, it will be outputted to several files
			if(type) for(int typecheck=0; typecheck < par_noEA; typecheck++) if(type & ( 1 << typecheck ) ) {
				//output replicators
				files[typecheck] << seq 
					<< '\t' << replicator.get_str() 
					<< '\t'	<< replicator.get_mfe() 
					<< '\t'	<< replicator.getPfold() 
					<< '\t'	<< replicator.Pdeg 
					<< '\t'	<< replicator.get_no_sites() 
					<< '\t'	<< replicator.getR() 
					<< '\t' << seq.length() //note: in contrary of usual output, here I output length instead of M (M has no meaning here) 
					<< '\t'	<< replicator.get_type(); 

				for(double *a = replicator.geta(), *a_until = replicator.geta() + par_noEA; a != a_until; a++){
					files[typecheck] << '\t' << *a;
				}

				files[typecheck] << std::endl; 			

				//output complement replicators
				compl_replicator.replicate( replicator );
				compl_files[typecheck] << *(compl_replicator.get_seq())  
					<< '\t' << compl_replicator.get_str() 
					<< '\t'	<< compl_replicator.get_mfe() 
					<< '\t'	<< compl_replicator.getPfold() 
					<< '\t'	<< compl_replicator.Pdeg 
					<< '\t'	<< compl_replicator.get_no_sites() 
					<< '\t'	<< compl_replicator.getR() 
					<< '\t' << compl_replicator.get_seq()->length() //note: in contrary of usual output, here I output length instead of M (M has no meaning here) 
					<< '\t'	<< compl_replicator.get_type(); 

				for(double *a = compl_replicator.geta(), *a_until = compl_replicator.geta() + par_noEA; a != a_until; a++){
					compl_files[typecheck] << '\t' << *a;
				}

				compl_files[typecheck] << std::endl; 			
			} //output complementer strand too
//		}
//		else {
//			std::cerr << "ERROR: calculation of number of types (" << type << " - max: " << no_types << ") were wrong!" << std::endl;
//		}
		
		//die
		replicator.die();
		compl_replicator.die();
	}








	for(auto f = files.begin(); f != files.end(); f++) (*f).close();
	for(auto f = compl_files.begin(); f != compl_files.end(); f++) (*f).close();
	//close rng
	gsl_rng_free(r);
	return 0;
}


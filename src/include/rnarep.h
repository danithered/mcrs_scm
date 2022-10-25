#ifndef _RNAREP_
#define _RNAREP_

#include <string.h>
//#include <iostream>
#include <vector>
#include <cmath>

extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
//#include <ViennaRNA/fold_vars.h>
}

#include "randomgen.h"
#include "parameters.h"
#include "dv_tools.h"
#include "annot.h"

namespace rnarep {
	//values
	extern char bases[6];

	//functions
	double length_depCalc(int n);
	double m_sigmaCalc(int m);
	
	char RNAc2cc(char rna);
	int RNAc2i(char rnachar);

	//classes
	class CellContent{
		public:
			//double value1; //temporary, just for testing

			double Pdeg; //degradation rate
			bool empty; // is this cell empty or not


			static int no_replicators; //number of replicators

			static dvtools::quickPosVals length_dep;
			static dvtools::quickPosVals m_sigma;
			static dv_annot::PatternPool patterns;

			//initialiser
			CellContent():Pdeg(0), empty(true), type(0), annot_level(0){	
				// allocate memory for MFE structure (length + 1)
				//str = (char *) vrna_alloc(sizeof(char) * ( MAXLEN  + 1));
				str = new char [MAXLEN + 1];

				// allocate memory for enzymaitc activities
				a = new double [par_noEA];

				seq.reserve(MAXLEN);
				//ins.reserve(MAXLEN+1);
				//subs.reserve(MAXLEN+1);
				//dels.reserve(MAXLEN);


				int len=0;
				
				//initialise array for seq
				//if ( !(seq = new char [MAXLEN]) ) std::cerr << "ERROR: could not allocate space for seq" << std::endl;
				//seq = new char [MAXLEN]; 
				//std::memset(seq, '\0', MAXLEN);

				//add random values for seq
				if(gsl_rng_uniform(r) < par_init_grid) {
					len = gsl_rng_uniform_int(r, 10) + 1; //initiate with n bases, n is [1,10] equal distribution
					for(int b=0; b < len; b++){
						seq.push_back( rnarep::bases[gsl_rng_uniform_int(r,4)] );
					}
					annotate();
				}
				else { //empty
					//empty = false; //to die properly
					die();
				}


				
//				std::cout << "initialised without content" << par_noEA << std::endl;
			 }

			CellContent(std::string input_str){
				empty=true;
				seq.reserve(MAXLEN);
				//ins.reserve(MAXLEN+1);
				//subs.reserve(MAXLEN+1);
				//dels.reserve(MAXLEN);

				seq = input_str;
				// allocate memory for MFE structure (length + 1)

				//str = (char *) vrna_alloc(sizeof(char) * ( MAXLEN  + 1));
				str = new char [MAXLEN + 1];

				// allocate memory for enzymaitc activities
				a = new double [par_noEA];

//				std::cout << "initialised" << par_noEA << std::endl;
	
				if(seq.length()){
					annotate();
				}
				else {
					die();
				}
			}

			~CellContent(){
				//delete [] (seq);
				//cleanup memory
				delete [] (str);
				delete [] (a);
			}

			void operator =( std::string& templ){
//				std::cout << "CellContent assignment from templ " << templ << " on empty (" << empty << ") cell.\n";				
				die(); // !empty is checked in die()

				//check if new seq is ok
				if(templ.length() && (templ != "N") && (templ != "0") ){ //length is not zero AND seq is not "N" or "0"
					seq = templ;
					//if(seq.length()) 
					annotate();
//					std::cout << "added new replicator " << seq << std::endl;
				}
//				else std::cout << "cell stayed empty" << std::endl;
			}

			void die();

			void replicate( CellContent &templ);
			
			void replicate_clear( CellContent &templ);

			double geta(int no);

			double* geta();

			double getR();

			int get_no_sites();

			int get_no_acts();

			unsigned long long int get_type();

			unsigned long long int get_type_rev();

			unsigned long long int get_prev_type();

			int get_length();

			double get_mfe();

			std::string * get_seq();

			char * get_str();

			double getPfold();

		//private:
			char *str;
			float mfe;
			double Pfold;
			int no_sites;
			int no_acts;
			std::string seq; //sequence of replicator
			double R; // replication rate
			double *a; //enzymatic activities
			unsigned long long int type; //an integer indicating the enzimatic activities that can be found in a replicator
			unsigned long long int prev_type; //an integer indicating the enzimatic activities that can be found in a replicator
			//std::vector<int> ins;
			//std::vector<int> subs;
			//std::vector<int> dels;

			int annot_level;

			//FUNCTIONS
			void annotate() { //getting mfe and Pdeg. Have to be called every time the sequence changes (and not for empty)
				empty = false;
				annot_level = 1;
				prev_type = 0;
				
				// predict Minmum Free Energy and corresponding secondary structure
				//mfe = vrna_fold(seq.c_str(), str);
				mfe = (double) vrna_fold(seq.c_str(), str);
					  
				//calculate Pdeg
				//Pdeg = 0.9 - 0.8 * (mfe<par_Emin?par_Emin:mfe) / par_Emin ;
				//Pdeg = par_maxPdeg -  par_rangePdeg* (mfe<par_Emin?par_Emin:mfe) / par_Emin ;
				//Pdeg = 0.05 - 0.049 * ( (mfe<par_Emin)?par_Emin:mfe) / par_Emin ;
				//Pdeg=0.07;
				Pdeg = par_minPdeg + par_rangePdeg * std::exp(par_flexPdeg * mfe);
				
				no_replicators++;
			} 

			void annotate2() { //getting Pfold, activities and type
				annot_level = 2;

				type = 0;
				no_acts= 0;

				//calculate Pfold
				// exp(-cE) / (1 + exp(-ce)) = 1 - 1 / (1 + exp(-cE)) , E = MFE, c = parameter
				Pfold = 1 - 1/( 1 + std::exp( par_c * mfe) );
				
				//annotata
			 	no_sites = patterns.search((char*) seq.c_str(), str, a);

//				if(no_sites) std::cout << "pattern.search has found an enzyme" << std::endl; else std::cout << "pattern.search has found a parazyte" << std::endl;

				//compute a from alpha
				for(int act = 0; act < par_noEA; act++) {
					if(a[act]){ //if there is such an ezymatic activity
						//unsigned long long int one;
						//one = 1;
						a[act] = Pfold * a[act] * m_sigma[no_sites]; //mind that m_sigma in now reciproc of prev function! 
						//adding to type
						type += ( ( (unsigned long long int) 1) << act);
						no_acts++; //note, that in this version return value of search is number of motifs, not number of activities!!
					}
					else a[act] = 0.0; // if there is no such activity
				}

			}

			void annotate3(){ //getting R
				annot_level = 3;
				//calculate R
				//R = g / (b1 + b2*L) * (l + (1 - Pfold)) , where L and Pfold are variables
				//if W(L) = g / (b1 + b2*L) , and W[L] is precalculated (L is int)
				//and ll = l + 1
				//R = W(L) * ( ll - Pfold )
				R = length_dep[ seq.length() ] * (par_ll - Pfold) ;

			}

	};	
}


#endif


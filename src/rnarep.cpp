#include "rnarep.h" 


namespace rnarep {
	char bases[6] = "GCUAN";

	char RNAc2cc(char rna){
		/* Function by: Enrico Sandro Colizzi */
		switch (rna) {
			case 'A':
				return 'U';
			case 'C':
			  	return 'G';
			case 'G':
				return 'C';
			case 'U':
				return 'A';
		      	default:
			        //fprintf(stderr,"RNAc2cc: I got a non RNA character (%c)\n", rna);
				std::cerr << "RNAc2cc: I got a non RNA character (" << rna <<")" << std::endl;
				return '\0';
		}
			
	}

	int RNAc2i(char rnachar){
		switch (rnachar) {
			case 'G':
			        return 0;
			case 'C':
				return 1;
			case 'U':
				return 2;
			case 'A':
				return 3;
			default:
				std::cerr << "RNAc2i: I got a non RNA character (" << rnachar << ")" << std::endl;
				return -1;
		}
	}


	double length_depCalc(int L){
		return par_g / (par_b1 + par_b2 * L);
	}

	double m_sigmaCalc(int m){
		return 1/(std::pow( (double) m, par_sigma));
	}

	dvtools::quickPosVals CellContent::length_dep(100, &length_depCalc);
	dvtools::quickPosVals CellContent::m_sigma(par_noEA*3, &m_sigmaCalc);

	int CellContent::no_replicators;

	dv_annot::PatternPool CellContent::patterns;

	//Functions for CellContent
	
	void CellContent::die(){
		if (!empty) {
			empty = true;
			type = 0;
			seq.clear();
			Pdeg = 0;
			//R = 0;
			annot_level = 0;
			no_replicators--;
			//for (int act = 0; act < par_noEA; act++) a[act] = 0;
		}
	}

	void CellContent::replicate_clear(CellContent &templ){
		if(!empty) die();
		for(auto old_it = templ.seq.rbegin(); old_it != templ.seq.rend(); old_it++){
				seq.push_back( RNAc2cc( (char) *old_it) ); //good (correct copying)
		}

		//for now it is seq is added -> need to annotate!!
		if(seq.length()){
			annotate();
			prev_type = templ.get_type(); 
		}
		/*else {
			die();
		}*/
	}

	void CellContent::replicate(CellContent &templ){
		
		for(auto old_it = templ.seq.rbegin(); old_it != templ.seq.rend(); old_it++){
			if( gsl_rng_uniform(r) < par_insertion ) seq.push_back( bases[gsl_rng_uniform_int(r, 4)] ); //if there is an insertion add random base
			if( gsl_rng_uniform(r) > par_deletion ) { //if there is no deletion copy template
				if( gsl_rng_uniform(r) < par_substitution ) seq.push_back( bases[( RNAc2i( RNAc2cc( (char) *old_it) ) + gsl_rng_uniform_int(r, 3) + 1) % 4] ); //wrong (substitution)
				else seq.push_back( RNAc2cc( (char) *old_it) ); //good (correct copying)
			}
		}
		if( gsl_rng_uniform(r) < par_insertion ) seq.push_back( bases[gsl_rng_uniform_int(r, 4)] ); //if there is an insertion add random base to its end
		


		//for now it is seq is added -> need to annotate!!
		if(seq.length()){
			annotate();
			prev_type = templ.get_type(); 
			/* for testing purposes: dont forget to comment it out!! */
//			if(prev_type == 0 && type == 0) die(); //it will cause perasites to die immediately - WARNING: just for testing purposes, comment it out for main simulations!!!"
			/*********************************************************/
		}
		else {
			die();
		}
		

	}

	double CellContent::geta(int no){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return 0;
		}
		
		return a[no];
	}

	double* CellContent::geta(){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return NULL;
		}
		
		return a;
	}

	double CellContent::getR(){
		if(annot_level < 3) { 
			if (annot_level < 2){
				if (annot_level) annotate2();
				else return -1.0 ;
			}
			annotate3();
		}

		return R;

	}

	int CellContent::get_no_sites(){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return -1;
		}
		
		return no_sites;
	}

	int CellContent::get_no_acts(){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return -1;
		}
		
		return no_acts;
	}

	unsigned long long int CellContent::get_type(){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return -1;
		}
		
		return type;

	}

	unsigned long long int CellContent::get_type_rev(){
		CellContent cell;
		cell.replicate_clear(*this);
		
		unsigned long long int output = cell.get_type();
		cell.die();

		return output;
	}

	unsigned long long int CellContent::get_prev_type(){
		if (annot_level < 2){ 
			if (annot_level) annotate2();
			else return -1;
		}
		
		return prev_type;

	}

	int CellContent::get_length(){
		if (empty) return 0;
		else return seq.length();
	}

	double CellContent::get_mfe(){
		if(annot_level < 1) return 0;
		else return mfe;
	}

	std::string * CellContent::get_seq(){
		if(empty) return NULL;
		else return &seq;
	}

	char * CellContent::get_str(){
		if(empty) return NULL;
		else return str;
	}

	double CellContent::getPfold(){
		if (annot_level < 2){
			if (annot_level) annotate2();
			else return -1;
		}
		
		return Pfold;

	}
}


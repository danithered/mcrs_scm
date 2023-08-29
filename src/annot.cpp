#include "annot.h"


namespace dv_annot{


	dvtools::quickPosVals PatternPool::gcBonusFunc(17, [](int x){return 0.5 + 1/(1+std::exp(par_gc_bonus * (double)x)); });

	//Functions for Subrule
	void Subrule::addBases(int _no_bases, char *_base, int *_pos) {
		no_bases = _no_bases;
		base = new char[no_bases];
		pos = new int[no_bases];
		std::strncpy(base, _base, no_bases);
		std::memcpy(pos, _pos, sizeof(int) * no_bases);

		//count GCs in subrule
		no_GC_in_pattern = 0;
		for(int i = 0; i < no_bases; i++){
			if(base[i] == 'G' || base[i] == 'C') no_GC_in_pattern++;
		}
	}

	//Funcion to initialise base
	//initBases(int _no_bases) {
	//	no_bases = _no_bases;
	//	base = new char[no_bases];
	//	pos = new int[no_bases];
	//}

	//Function to add base
	//void setBase (int number, int pos, char base){
	//	pos[number] = pos;
	//	base[number] = base;
	//}
	

	//Functions for Rule
	void Rule::setPattern(std::string p){
		pattern_length = p.length();
		if(pattern_length){
			pattern = new char [pattern_length + 1];
			std::strcpy(pattern , p.c_str());
			//pattern[pattern_length] = '\0';
		}
//				std::cout << "Rule pattern set" << std::endl;
	}

	void Rule::setSubrules(std::ifstream &file){
//				std::cout << file.rdbuf();
		std::string line;
		std::vector<char> chars;
		std::vector<int> positions;
		char c;
		int p = -1;

		for( int sr=0; sr < no_subrules; sr++){
//					if(! (file) ) std::cout << "File is NULL at sr= " << sr << std::endl;
			std::getline(file, line); //jump to next line
			std::getline(file, line);
			std::istringstream linestream(line);

//					std::cout << "Line read: " << line << std::endl;
			chars.clear();
			positions.clear();

			while(linestream >> p) { //read in position
				linestream >> c; //read in character
				chars.push_back(c);
				positions.push_back(p - 1); //indexing i c starts from 0, while input file starts from 1
//						std::cout << "new subrule element added: " << c << "at pos " << p << std::endl;
			}
			subrules[sr].addBases( chars.size(), chars.data(), positions.data()); //.data() is only cpp11, replace with &x[0] 
			//adding values
//					if(! (file) ) std::cout << "File is NULL at (2) sr= " << sr << std::endl;
			for(int eaval = 0; eaval < par_noEA; eaval++) {
//						if(! (file) ) std::cout << "File is NULL at (3, eaval= " << eaval << ") sr= " << sr << std::endl;
				file >> subrules[sr].value[eaval];
//						std::cout << "value read (at eaval" << eaval << "):" << subrules[sr].value[eaval] << std::endl;
			}
//					if(! (file) ) std::cout << "File is NULL at (4) sr= " << sr << std::endl;
		}
//				std::cout << "Rule: subrule added" << std::endl;
	}

	//Functions for PatternPool
	void PatternPool::readFile(char *filename, const char *copy){
		PatternPool::readFile(filename);
		if(std::strlen(copy)){
			std::filesystem::copy(filename, copy);
		}
	}
	void PatternPool::readFile(char *filename){
		int no = 0; 
		int min = 0; //number of ea, subrules
		std::string word;
		
		std::ifstream file;
		file.open(filename);
		if(!file.is_open()) std::cerr << "ERROR: file (" << filename << ") not found!" << std::endl;
//				std::cout << file.rdbuf();

		//read first word which is number of enzymatic activities
		file >> no;
		if( no != par_noEA){
			//if (par_noEA) delete [] (a);
			par_noEA = no;
			//if (par_noEA) a = new double [par_noEA];
			//for(int i=0; i < par_noEA; i++) a[i] = 0.0;
		}
		no = 0;
		no_sites = new int[par_noEA];

		//read in rest of file
		while (file >> no){ //read number of subrules
			rules.push_back(Rule(no)); //created Rule
			file >> word; //read pattern
//					std::cout << "word read: " << word << std::endl;
			rules.back().setPattern(word);
			rules.back().setSubrules(file); //set subrules (all)
		}
//				if(file )std::cout << file.tellg() << " " << file.gcount() << std::endl; else std::cout << "file is NULL" << std::endl; 

		//closing file
		file.close();

		//sorting rules by pattern_length 
		if(rules.size() > 1) {
			class Rule memory(0);
			for( no = 0; (unsigned int) no < (rules.size() - 1); no++){
				for(unsigned int order = no + 1, min = no; order < rules.size(); order++){
					if(rules[order].pattern_length < rules[min].pattern_length ) min = order;
				}
				if(min != no) {
//							std::cout << sizeof(memory) << " " << sizeof(rules[min]) << " " << sizeof(rules[no]) << std::endl;
					/*
					memcpy(&memory, &rules[min], sizeof(rules[min]) ); //save min to memory
					memcpy(&rules[min], &rules[no], sizeof(rules[no]) ); //copy no to min
					memcpy(&rules[no], &memory, sizeof(memory) ); //copy memory to no 
					*/
					memory = rules[min]; //save min to memory
					rules[min] = rules[no]; //copy no to min
					rules[no] = memory; //copy memory to no
				}

			}
			memory.pattern_length = 0; //to ensure no double free will happen
			memory.no_subrules = 0; //detto
		}

	}
	
	//void clear(){
	//	//if(a) for(int i=0; i < par_noEA; i++) a[i] = 0.0;
	//	sites = 0;
	//}

	int PatternPool::search(char *seq, char *str, double *acts){
		char *templ, *templ_seq;
		int templ_length = std::strlen(seq);
		double stepwise[] = {0.0, 0.1, 0.8, 1.0};

		//clear();
		int sites=0;
		for(int act = 0; act < par_noEA; act++) {
			no_sites[act] = 0;
			acts[act] = 0;
		}
		
//		std::cout << "search this: " << seq << "\t" << str << "\t a:"; for (int pr = 0; pr < par_noEA; pr++ ) std::cout << acts[pr] << " "; std::cout << std::endl;
		
		if(templ_length) for(unsigned int search = 0; search < rules.size() && rules[search].pattern_length <= templ_length ; search++){ // goes tru rules
//			std::cout << "search = " << search << std::endl;

			for(templ = std::strstr(str, rules[search].pattern) ; templ != NULL; templ = std::strstr(++templ, rules[search].pattern)){ // finds rule's patterns
//				std::cout << "looking for structure [" << rules[search].pattern << "] in str (search = " << search << ")"  << std::endl;
				//templ is a pointer to the start of the str pattern
				//look for the subrules
				for(int sr = 0, sr_max = rules[search].no_subrules; sr < sr_max; sr++){ // if pattern found goes tru subrules
//					std::cout << "sr = " << sr << std::endl;
					//looking tru the subrules
					char *base = rules[search].subrules[sr].base;
					int *pos = rules[search].subrules[sr].pos;
					int *pos_max = pos + rules[search].subrules[sr].no_bases;
					double subrules_apply = 0.0;

//					std::cout << "compare " << *(seq + (templ - str + *pos)) << "and" << *base << " (from " << rules[search].subrules[sr].no_bases << " bases)" << std::endl;
					//for(templ_seq = seq + (templ - str); pos != pos_max && templ_seq[*pos] == *base; base++, pos++){} //checks subrule
					for(templ_seq = seq + (templ - str); pos != pos_max; base++, pos++){
						if(templ_seq[*pos] == *base) subrules_apply++; //checks subrule
					}
//					std::cout << "pos " << pos << " pos max "<< pos_max << std::endl;
//
					//if number of subrules = 3 -> go stepwise like balazs, else go linearly
					if(rules[search].no_bases == 3){
						subrules_apply = stepwise[(int) subrules_apply];
					} else {
						subrules_apply /= (double) rules[search].no_subrules;
					}

					if (subrules_apply > 0){ //subrule applies
						sites++;

						//looking tru activities in subrule
						for(int act = 0; act < par_noEA; act++){
//							std::cout << rules[search].subrules[sr].value[act] << " added to activity " << act << std::endl;
							double curr_act = rules[search].subrules[sr].value[act];

							if(curr_act){ //if subrule adds a valid activity
								/******************************/
								// now find the uc ratio!
								if(par_gc_bonus) {
									//int open_bases = 0; //number of points
									int gc_num = 0; //it will compute the ratio of pirimidin bases in open "." positions in the structure

									//compute number of open bases in motif and number of C/Us in open places
									for(int b = 0; b < rules[search].pattern_length; b++){
										if(templ[b] == '.' && (templ_seq[b] == 'G' || templ_seq[b] == 'C')){
											//open_bases++;
											//if(templ_seq[b] == 'G' || templ_seq[b] == 'C') 
											gc_num++;
										}
									}

									//calc
									gc_num -= rules[search].subrules[sr].no_GC_in_pattern; //substract GC that is there as a part of the subrule
									if( gc_num > 0) {
										//acts[act] += curr_act * (1 + par_gc_bonus * (double) gc_num / open_bases); // if there is bonus system and bonus
										acts[act] += subrules_apply * curr_act * gcBonusFunc[gc_num]; // if there is bonus system and bonus
										 
									} else {
										acts[act] += subrules_apply * curr_act; // if there is bonus system, but no bonus (only core activity)
									}
								} else {
									//no bonus system
									acts[act] += subrules_apply * curr_act;
								}

								/******************************/
								no_sites[act]++;
							}
						}
						break; // if subrule applies, doesnt checks for other subrules in this pattern -> goes to find next pattern
					}
					// it is the wanted structure

				} //sr in subrules
			} // pattern search in sequence
		} // search in rules

		/*//calculating average activity values and number of activities
		for(int act = 0; act < par_noEA; act++){
			if(no_sites[act]){
				sites++;
				acts[act] = acts[act] / no_sites[act]; // average of all the same activities (DV strongly disagrees with this)
			}
		}*/

//		std::cout << "search this: " << seq << "\t" << str << "\t a:"; for (int pr = 0; pr < par_noEA; pr++ ) std::cout << acts[pr] << " "; std::cout << std::endl;

		return(sites); //currently it is number of motifs in the structure!! It is used later too, mind if you change it!
	}

	void PatternPool::printRules(){
		for(unsigned int r = 0; r < rules.size(); r++){
			std::cout << "############################" << std::endl;
			std::cout << "RULE no " << r+1 << std::endl << std::endl;
			std::cout << "pattern: " <<  rules[r].pattern << std::endl << std::endl;
			std::cout << "SUBRULES:" << std::endl;
			for(int sr = 0; sr < rules[r].no_subrules; sr++){
				std::cout << "\t" << "Subrule no " << sr+1 << " from " << rules[r].no_subrules << std::endl;

				std::cout << "\t";
				for(int b =0; b < rules[r].subrules[sr].no_bases; b++){
					std::cout << "\t" << rules[r].subrules[sr].pos[b]+1 << " " << rules[r].subrules[sr].base[b];
				}
				std::cout << std::endl << "\tvalues: ";
				for(int v = 0; v < par_noEA; v++){
					std::cout << rules[r].subrules[sr].value[v] << " ";
				}
				std::cout << std::endl << std::endl;
			}
		}
	}

}


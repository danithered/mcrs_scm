#ifndef _SEQANNOT_
#define _SEQANNOT_

#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <filesystem>

#include "parameters.h"
#include "dv_tools.h"

namespace dv_annot{

	class Subrule {
		public:
			double *value; //the value of the structure
			int no_bases; //number of bases as subrule
			int no_GC_in_pattern;

			char *base;
			int *pos;

			//Constructor
			Subrule(){
				no_bases = 0;
				no_GC_in_pattern = 0;
				value = new double [par_noEA];
//				std::cout << "Subrule allocated (empty constr) " << par_noEA << " values" << std::endl;
//				std::cout << "Subrule initialised" << std::endl;
			}

			//Copy Constructor
			Subrule(const Subrule &obj){
				int i = 0;
				no_bases = obj.no_bases;
				no_GC_in_pattern = obj.no_GC_in_pattern;
				if(no_bases){
					base = new char[no_bases];
					pos = new int[no_bases];
					for(i = 0; i < no_bases; i++){
						base[i] = obj.base[i];
						pos[i] = obj.pos[i];
					}
				}
				value = new double [par_noEA];
//				std::cout << "Subrule allocated (copy constr) " << par_noEA << " values" << std::endl;
				for(i = 0; i < par_noEA; i++){
					value[i] = obj.value[i];
				}
//				std::cout << "Subrule copied" << std::endl;
			}

			//Operator =
			void operator = (const Subrule &obj){
				int i = 0;
				no_bases = obj.no_bases;
				no_GC_in_pattern = obj.no_GC_in_pattern;
				if(no_bases){
					base = new char[no_bases];
					pos = new int[no_bases];
					for(i = 0; i < no_bases; i++){
						base[i] = obj.base[i];
						pos[i] = obj.pos[i];
					}
				}
				//reallocate memory - probably completely unnecessary, but it cost nothing
				delete [] (value);
				value = new double [par_noEA];
//				std::cout << "Subrule allocated (operator =) " << par_noEA << " values" << std::endl;

				//copy value
				for(i = 0; i < par_noEA; i++){
					value[i] = obj.value[i];
				}
//				std::cout << "Subrule copied by = operator" << std::endl;
			}

			//Destructor
			~Subrule(){
//				std::cout << "Subrule destructor called" << std::endl;
				delete [] (value);
//				std::cout << "Subrule deleted " << par_noEA << " values" << std::endl;

				if (no_bases) {
					delete [] (base);
					delete [] (pos);
				}
//				std::cout << "Subrule destructor ended" << std::endl;
			}


			//Function to add all the bases at once
			void addBases(int _no_bases, char *_base, int *_pos);


	}; //Subrule

	class Rule {
		public:
			int pattern_length;
			char *pattern;

			int no_subrules;
			class Subrule *subrules;

			//Constructor
			Rule(int _no_subrules) : no_subrules{_no_subrules} {
//				std::cout << "Rule init started" << std::endl;
				if(no_subrules) subrules = new class Subrule [no_subrules];
				pattern_length = 0;
//				std::cout << "Initialised Rule" << std::endl;
			}

			//Copy Constructor
			Rule(const Rule &obj){
//				std::cout << "Rule copy init started" << std::endl;
				pattern_length = obj.pattern_length;
				if(pattern_length){
					pattern = new char [pattern_length + 1];
					std::strcpy(pattern , obj.pattern);
				}

				no_subrules = obj.no_subrules;
				if(no_subrules) {
					subrules = new class Subrule [no_subrules];
					for (int s = 0; s < no_subrules; s++){
						subrules[s] = obj.subrules[s]; 
					}
				
				}
				else {
					subrules = NULL;
				}
//				std::cout << "Rule Copied" << std::endl;
			}

			//Destructor
			~Rule(){
//				std::cout << "Rule (" << no_subrules << ") destructor called" << std::endl;
				if(no_subrules) delete [] (subrules);
//				std::cout << "Rule destructor halfway" << std::endl;
				if (pattern_length) delete [] (pattern);
//				std::cout << "Rule destructor ended" << std::endl;
			}


			//Functions
			void setPattern(std::string p);

			void setSubrules(std::ifstream &file);

	}; //Rule

	class PatternPool {
		public:
			//extern int par_noEA;
			std::vector<class Rule> rules;	//pointer for rules
			//double *a;		//enzim activities for a search
			//int sites;

			//Constructor
			PatternPool(){
				//no_sites = new int[par_noEA];
//				std::cout << "Pattern pool init started" << std::endl;
				//sites = 0;
				//if (par_noEA) a = new double [par_noEA];
				//for(int i=0; i < par_noEA; i++) a[i] = 0.0;
				//std::memset(a, 0, sizeof(double) * par_noEA); //remember: use memset only for 0 in case of int!!
//				std::cout << "Pattern pool init ended" << std::endl;
			}

			//Destructor
			~PatternPool(){
//				std::cout << "PatternPool destructor called" << std::endl;
				//if (par_noEA) delete [] (a);
				rules.clear();
				delete [] (no_sites);
//				std::cout << "PatternPool destructor ended" << std::endl;
			}

			//file input
			void readFile(char *filename);
			void readFile(char *filename, const char *copy);
			
			//clearing activities
			//void clear();

			//searching pattern
			int search(char *seq, char *str, double *acts);

			//print out rules   DONT FORGET TO COMMENT IT OUT!!!
			void printRules();

		private:
			int *no_sites;
			static dvtools::quickPosVals gcBonusFunc;



	}; //PatternPool

}

#endif


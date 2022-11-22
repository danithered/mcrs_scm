#ifndef _CADV_
#define _CADV_

#include <vector>
#include "randomgen.h"
#include "dv_tools.h"
#include "rnarep.h"
#include <cmath>
#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include <list>

#define SINGLESQ 2.0
#define SINGLEHEX 3.0
#define MOORE_NEIGH 4.0
#define VON_NEUMANN_NEIGH 3.0
#define MARGOLUS_NEIGH -1.0
#define NEIGH5X5 8.0
#define HEX1 4.0


namespace cmdv {
	class Compart{
		public:
			//Compart content
			std::list<class rnarep::CellContent> reps;

			class CompartPool *parent;
			//bool updateable;
			bool alive;
			static unsigned int no_alive;

			//Functions
			//Compart base functions
			Compart(){
				parent = NULL;
				alive=false;
				//metabolism = 0;
				leftover = 0;
				//updateable = true;
				reciproc_noEA = 1 / (double) par_noEA;
			}

			~Compart(){
				//for(auto repdel = reps.begin(); repdel != reps.end(); repdel++) delete *repdel;
				//for(auto repdel = wastebin.begin(); repdel != wastebin.end(); repdel++) delete *repdel;
				
			}

			///owerwrite one compart with other
			void operator =(Compart& origin){
				clear();
				//for(auto rep = origin.reps.begin(); rep != origin.reps.end(); rep++){
					reps.assign( origin.reps.begin(), origin.reps.end() ); 
				//}
				//no need to delete temp_compart, split() will do that anyway...
			}

			//add replicator
			std::list<rnarep::CellContent>::iterator add(rnarep::CellContent rep);
			std::list<rnarep::CellContent>::iterator add();
			std::list<rnarep::CellContent>::iterator add(std::string newseq);
			std::list<rnarep::CellContent>::iterator add(std::list<rnarep::CellContent>::iterator it, std::list<rnarep::CellContent> &from);

			//kill replicator (CellContents own die() and storing it in wastbin)
			void die(std::list<rnarep::CellContent>::iterator rep);

			void sleep();
			void wake();

			//calculate metabolism around replicator
			double M();

			//an update step on this cell
			void update();
			std::list<rnarep::CellContent>::iterator replication();

			//split to two compartments
			int split();

			//clear to ensure it can be rewritten during Moran process (other split)
			void clear();


		//private:
			static std::list<class rnarep::CellContent> wastebin;
			double reciproc_noEA;
			double leftover;
	
	};

	
	class CompartPool {
		public:
			unsigned int size;
			int time;
			//int no_replicators;

			Compart **comparts;
			Compart **temp_comparts;
			unsigned int used_temp;
			
			std::string savedir;

			//OUTPUTS
			std::ofstream output;

			//FUNCTIONS

			//Constructor 1
			CompartPool(int _size=300): size(_size), used_temp(0){
				time=0;
				//saving_freq = 0;

				comparts = new class Compart* [size];
				for(Compart **comp = comparts, **endcomp = comparts+size; comp != endcomp; comp++){
					*comp = new class Compart;
					(*comp)->parent = this;
				}

				//maximum number of new comparts is equal to number of comparts
				temp_comparts = new class Compart* [size*10];
				for(Compart **comp = temp_comparts, **endcomp = temp_comparts+size*10; comp != endcomp; comp++){
					*comp = new class Compart;
					(*comp)->parent = this;
				}

				//no_replicators = 0;

//				std::cout << "Basic Constructor Called" << std::endl;
			}
			
			//Constructor 2 - for deserialisation
			CompartPool() {}
			
			//Destructor
			~CompartPool(){
				if(size) {
					for(unsigned int i = 0; i < size; ++i ){
						delete comparts[i];
						delete temp_comparts[i];
					}
					delete [] (comparts);
					delete [] (temp_comparts);
				}
//				std::cout << "Deconstructor Called" << std::endl;
			}
			
			//virtual void Update(int cell);
			
			//Operators
			
			
			//Functions
			
			///gives back pointer to nth comp
			inline Compart* get(int cell) {
				return( *(comparts + cell) );
			}

			///gives back pointer to random comp
			inline Compart** get() {
				return( comparts + gsl_rng_uniform_int(r, size) );
			}
			
			///gives back pointer to random comp (excluded arg)
			inline Compart* get(Compart *except) {
				if(size < 2) return NULL;

				int n = (int) (except - *comparts);
				return( *(comparts + (n + gsl_rng_uniform_int(r, size-1) + 1) % size) );

				/*compart* out = comparts + gsl_rng_uniform_int(r, size);
				while(except == out) comparts + gsl_rng_uniform_int(r, size);
				return(out);*/
			}
			
			///initialises matrix with predefined values, randomly
			//void init(std::string* pool, double* probs, int no_choices); 
			
			///initialises matrix from textfile
			void init_fromfile(char * infile); 

			//bool compartFromFile(const char * infile);
			
			//clear updateable flag
			//void all_updateable();
			
			//Updates

			///One update step - DOES NOT WORK
			//int updateStep(int cell);
			
			// Throwing back comparts from temp
			void compartShower();

			///Random update
			int rUpdate(int gens);

			//Update according to a random order (in every generation all cells will be updated)
			int oUpdate(int gens);

			// Outputs
			// return values: 
			// 0 (OK), 
			// 1 (could not find a new directory, most likely more than one simulations run with the same ID and seed), 
			// 2 (cant open output file)
			
			int openOutputs();

			void do_output();

			// return values: 0 (OK), 
			// 1 (no savedir specified), 
			// 2 (could not open file)
			int save();

		private:
			std::vector<int> out_no; //how many replicator has no act, act0, act1, etc.
			std::vector<int> out_noA; //how many replicator has alltogether 0, 1, 2, etc different activities
			std::vector<double> out_R; //mean R of replicators with no act, act0, act1, etc.
			std::vector<double> out_length; //mean length of replicators with no act, act0, act1, etc.
			std::vector<double> out_a; //mean activity of replicators with no act, act0, act1, etc. (the strength of the indicated activities of course)
			std::vector<double> out_mfe; //mean mfe of replicators with no act, act0, act1, etc.
	};
	
}

#include "cm_serialise.h"


#endif

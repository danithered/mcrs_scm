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
			std::vector<double> claims;

			class CompartPool *parent;

			//Functions
			//Compart base functions
			Compart(){
				parent = NULL;
				metabolism = 0;
			}

			~Compart(){
				//delete reps;
			}

			//add replicator
			int add();

			//kill replicator (CellContents own die() and storing it in wastbin)
			int die(rnarep::CellContent *rep);

			//calculate metabolism around replicator
			double M();

			//an update step on this cell
			void update();
			void replication();
			void degradation();

			//split to two compartments
			int split();

			//clear to ensure it can be rewritten during Moran process (other split)
			void clear();


		private:
			double metabolism; //have to invalidate after each update step 
			std::list<class rnarep::CellContent> reps;
			std::list<class rnarep::CellContent> wastebin;
	
	};

	
	class CompartPool {
		public:
			int size;
			int time;
			//int no_replicators;

			Compart *comparts;
			
			std::string savedir;

			//OUTPUTS
			std::ofstream output;

			//FUNCTIONS

			int grid_init() ;

			//Constructor 1
			CompartAut(int _size=300): size(_size){
				time=0;
				//saving_freq = 0;

				comparts = new class Compart [size];
				for(Compart *comp = comparts, *endcomp = comparts+size; compart != endcomp; compart++){
					compart->parent = this;
				}

				//no_replicators = 0;

//				std::cout << "Basic Constructor Called" << std::endl;
			}
			
			
			//Deconstructor
			~CompartPool(){
				if(size) delete [] (comparts);
//				std::cout << "Deconstructor Called" << std::endl;
			}
			
			//virtual void Update(int cell);
			
			//Functions
			
			///gives back pointer to nth comp
			inline Compart* get(int cell) {
				return(comparts + cell);
			}

			///gives back pointer to random comp
			inline Compart* get() {
				return(comparts + gsl_rng_uniform_int(r, size) );
			}
			
			///gives back pointer to random comp (excluded arg)
			inline Compart* get(Compart *except) {
				if(size < 2) return NULL;

				int n = (int) comparts - except;
				return( comparts + (n + gsl_rng_uniform_int(r, size-1) + 1) %% size );

				/*compart* out = comparts + gsl_rng_uniform_int(r, size);
				while(except == out) comparts + gsl_rng_uniform_int(r, size);
				return(out);*/
			}
			
			///initialises matrix with predefined values, randomly
			//void init(std::string* pool, double* probs, int no_choices); 
			
			///initialises matrix from textfile
			//void init_fromfile(char * infile); 
			
			//Updates

			///One update step - DOES NOT WORK
			int updateStep(int cell);

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


#endif

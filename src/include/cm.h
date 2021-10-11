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
			std::list<class rnarep::CellContent> reps;

			class CompartPool *parent;

			//Functions
			//Compart base functions
			Compart(){
				vals = new rnarep::CellContent; //creating place for cell content values

				no_met_neigh = no_repl_neigh = no_diff_neigh = 0; //no inic
			}

			~Compart(){
				if (no_repl_neigh) {
					delete [] (claims);
				}

				delete vals;
			}

			//calculate metabolism around replicator
			double M();

			//an update step on this cell
			void update();

			//give birth to another compartment
			int birth();

		private:
			double metablism;
	
	};

	
	class CompartPool {
		public:
			int size;
			int time;
			cadv::Ca_layout layout;
			
			Compart *matrix;
			
			std::string savedir;

			//OUTPUTS
			std::ofstream output;

			//FUNCTIONS

			int grid_init() ;

			//Constructor 1
			CompartAut(int size1=300, int size2=300, cadv::Ca_layout layout_type=square){
				time=0;
				diff = 0;
				//saving_freq = 0;
				nrow=size1;
				ncol=size2;
				layout = layout_type;

				rnarep::CellContent::no_replicators = 0;

				grid_init();

				if(!size) layout = empty;
				if(size1==1 || size2==1){
					if(size1==size2) {
						layout = single;
					}
					else {
						layout = line;
					}
				}

				
//				std::cout << "Basic Constructor Called" << std::endl;
			}
			
			//Constructor 2
			CompartPool(int size1, int size2, cadv::Ca_layout layout_type, Compart* pool, double* probs, int no_choices){
				CompartPool(size1, size2, layout_type);
				//init(pool, probs, no_choices);
			}
			
			//Deconstructor
			~CompartPool(){
				if(size) delete [] (matrix);
				//for(int i =0; i < neighbourhoods.size(); i++) {
				//	delete [] (neighbourhoods.at(i));
				//}
//				std::cout << "Deconstructor Called" << std::endl;
			}
			
			//virtual void Update(int cell);
			
			//Functions
			
			///gives back pointer to random cell
			inline Compart* rget(int cell) {
				return(matrix + cell);
			}
			
			///initialises matrix with predefined values, randomly
			void init(std::string* pool, double* probs, int no_choices); 
			
			///initialises matrix from textfile
			void init_fromfile(char * infile); 
			
			//Updates

			///One update step - DOES NOT WORK
			int updateStep(int cell);

			///Random update
			int rUpdate(int gens);

			//Update according to a random order (in every generation all cells will be updated)
			int oUpdate(int gens);

			// Outputs
			// return values: 0 (OK), 1 (could not find a new directory, most likely more than one simulations run with the same ID and seed), 2 (cant open output file)
			int openOutputs();

			void do_output();

			// return values: 0 (OK), 1 (no savedir specified), 2 (could not open file)
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

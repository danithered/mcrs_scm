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
			bool alive;
			bool awake; 
			static unsigned int no_alive;

			//Functions
			//Compart base functions
			Compart();

			~Compart(){}

			///owerwrite one compart with other
			void operator =(Compart& origin);

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
			std::vector<Compart *> temp_comparts;
			unsigned int used_temp;
			unsigned int no_last_splits;
			
			std::string savedir;

			//OUTPUTS
			std::ofstream output;

			//FUNCTIONS

			//Constructor 1
			CompartPool(int _size=300);
			
			//Constructor 2 - for deserialisation
			CompartPool() {}
			
			//Destructor
			~CompartPool();
			
			//Operators
			
			
			//Functions
			
			///gives back pointer to nth comp
			inline Compart* get(int cell);

			///gives back pointer to random comp
			inline Compart** get();
			
			///gives back pointer to random comp (excluded arg)
			inline Compart* get(Compart *except);
			
			///initialises matrix from textfile
			void init_fromfile(char * infile); 

			//bool compartFromFile(const char * infile);
			
			//Updates

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

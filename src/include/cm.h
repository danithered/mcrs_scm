#ifndef _CADV_
#define _CADV_

#include <vector>
#include "randomgen.h"
#include "dv_tools.h"
#include "rnarep.h"
#include "broken.hpp"
#include <cmath>
#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include <filesystem>
#include <list>
#include <map>
#include <unordered_set>

#define SINGLESQ 2.0
#define SINGLEHEX 3.0
#define MOORE_NEIGH 4.0
#define VON_NEUMANN_NEIGH 3.0
#define MARGOLUS_NEIGH -1.0
#define NEIGH5X5 8.0
#define HEX1 4.0

namespace fs = std::filesystem;

namespace cmdv {

	typedef std::map<unsigned int, std::list<std::string> > Bubbles;

	enum QuitCond{qnone=0, qreplicator=1, qcompart=2, qsplit=4, qfull=8, qalive=16};

	struct Outdata{
		int no; //how many replicator has no act, act0, act1, etc.
		double R; //mean R of replicators with no act, act0, act1, etc.
		double length; //mean length of replicators with no act, act0, act1, etc.
		double a; //mean activity of replicators with no act, act0, act1, etc. (the strength of the indicated activities of course)
		double mfe; //mean mfe of replicators with no act, act0, act1, etc.
		
		Outdata():no(0), R(0.0), length(0.0), a(0.0), mfe(0.0) {};

		void add(const double &_R, const double &_length, const double &_a, const double &_mfe ){
			++no;
			R += _R;
			length += _length;
			a += _a;
			mfe += _mfe;
		}
		
		double avg_R() const {return R/no;}
		double avg_length() const {return length/no;}
		double avg_a() const {return a/no;}
		double avg_mfe() const {return mfe/no;}
	};

	class Compart{
		public:
			class ScmRep : public rnarep::CellContent{
					Compart* vesicule;
					broken::FenwickNode<ScmRep*> *deathcount;
					broken::FenwickNode<ScmRep*> *repcount;
					
				public: 
					void setBindings(broken::FenwickNode<ScmRep*> *repBS, broken::FenwickNode<ScmRep*> *deathBS);
					Compart* assignCompart(Compart * comp);
					void updateDeg() const;
					void updateM() const;
					void updateRep(const double met) ;
					void degradets();
					void replicates();
					
					ScmRep(): rnarep::CellContent(), vesicule(nullptr), deathcount(nullptr), repcount(nullptr){}

			};

			//Compart content
			std::unordered_set<class ScmRep*> reps;

			class CompartPool *parent;
			static unsigned int no_alive;

			//Functions
			//Compart base functions

			//Constructor
			Compart();

			//Destructor
			~Compart(){}

			//Move constructor
			Compart(Compart && origin){std::cerr << "Undefined move constructor called" << std::endl;}
			///owerwrite one compart with other
			//void operator =(Compart& origin);

			//add replicator
			auto add(ScmRep* rep);
			ScmRep* add();
			ScmRep* add(std::string newseq);
			//std::list<ScmRep>::iterator add(std::list<ScmRep>::iterator it, std::list<ScmRep> &from);

			//kill replicator (CellContents own die() and storing it in wastbin)
			void die(ScmRep* rep);

			//calculate metabolism around replicator
			void refresh_M();
			void printReps();
			inline double get_M() const {return M;} 

			inline bool alive() {return _alive;}

			//an update step on this cell
			void replicate(Compart::ScmRep* const templ);

			//split to two compartments
			bool split();

			//clear to ensure it can be rewritten during Moran process (other split)
			void clear();

			double reciproc_noEA;
		private:
			double M;
			bool _alive;
	
	};

	
	class CompartPool {
		public:
			unsigned int size;
			int time;

			Compart **comparts;
			Compart::ScmRep *replicators;
			std::vector<Compart::ScmRep *> rep_stack;
			broken::reBrokenStick<Compart::ScmRep*> degpool;
			broken::reBrokenStick<Compart::ScmRep*> reppool;

			unsigned int no_last_splits; //< number of splits in last update step
			unsigned int no_last_replicates;
			unsigned int no_last_deaths;
			//unsigned int no_reps_last_in_alive; // to store the number of replicators in alive vesicules in previous generation
			
			std::string savedir;
			std::filesystem::path outpath;

			//OUTPUTS
			std::ofstream output;

			//FUNCTIONS

			//Constructor 1
			CompartPool(int _size=300);
			
			//Constructor 2 - for deserialisation
			CompartPool():degpool(1), reppool(1) {}
			
			//Constructor 3
			CompartPool(char *file);

			//Destructor
			~CompartPool();
			
			//Operators
			
			
			//Functions
			
			/// gives back pointer to nth comp
			inline Compart* get(int cell);

			/// gives back pointer to random comp
			inline Compart* get();
			
			/// gives back pointer to random comp (excluded arg)
			inline Compart* get(Compart *except);
			
			/// initialises matrix from textfile
			void init_fromfile(char * infile); 

			unsigned int discoverComparts(const char * sourcedir); 
			bool compartFromFile(const char * infile); // import compart to random place
			bool compartFromFile(const char * infile, const unsigned int n); // import compart to nth place
			void autoCompartInput();
			
			//Updates

			///Simple async update
			int Update(int gens);

			///Custom update
			int cUpdate(int gens);

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
			std::vector<Outdata> out_spec;
			std::vector<Outdata> out_gen;
			Outdata out_par;
			Outdata out_templ;
			std::vector<int> out_noA; //how many replicator has alltogether 1, 2, etc different activities (0 is par+templ)
						  
//			std::vector<int> out_no; //how many replicator has no act, act0, act1, etc.
//			std::vector<int> out_noA; //how many replicator has alltogether 0, 1, 2, etc different activities
//			std::vector<double> out_R; //mean R of replicators with no act, act0, act1, etc.
//			std::vector<double> out_length; //mean length of replicators with no act, act0, act1, etc.
//			std::vector<double> out_a; //mean activity of replicators with no act, act0, act1, etc. (the strength of the indicated activities of course)
//			std::vector<double> out_mfe; //mean mfe of replicators with no act, act0, act1, etc.

			Bubbles bubblefiles;
			unsigned int no_reps_last; // to store the number of replicators in previous generation

			bool testQuit() const; 
			unsigned int periodlength;
	};
	
}

#include "cm_serialise.h"


#endif

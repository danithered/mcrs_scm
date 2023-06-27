#include "cm.h"
#include "include/rnarep.h"

namespace cmdv {
	//int no_births=0;
	//int no_deaths=0;
	
	unsigned int Compart::no_alive = 0;
	
	inline void Compart::ScmRep::setBindings(broken::FenwickNode<ScmRep*> *repBS, broken::FenwickNode<ScmRep*> *deathBS){
		deathcount = deathBS;
		repcount = repBS;
	}

	// from: removes itself if it does not comes from stack
	// to: assotiates with it both ways
	inline Compart* Compart::ScmRep::assignCompart(Compart * comp){
		if(vesicule == comp) return vesicule;
		if(vesicule != nullptr){
			vesicule->reps.remove(this);
		}

		if(comp == nullptr){
			vesicule->parent->rep_stack.push_back(this);
		} else {
			comp->reps.push_back(this);
		}

		auto old_vesicule = vesicule;
		vesicule = comp;

		return(old_vesicule);

	}

	void Compart::ScmRep::updateDeg() const{
		deathcount->update(Pdeg);
	}

	void Compart::ScmRep::updateM() const{
		vesicule->refresh_M();
	}

	void Compart::ScmRep::updateRep(const double met) {
		repcount->update(met*getR());
	}

	Compart::Compart(): parent(NULL), reciproc_noEA(1 / (double) par_noEA), M(0.0), _alive(false){}

	/*void Compart::operator =(Compart& origin){
		clear();
		//for(auto rep = origin.reps.begin(); rep != origin.reps.end(); rep++){
			reps.assign( origin.reps.begin(), origin.reps.end() ); 
		//}
		//no need to delete temp_compart, split() will do that anyway...
	}*/

	void Compart::clear(){
		for(auto rep = reps.begin(), temp = rep; rep != reps.end(); temp=rep) {
			++rep;
			die(*temp);
		}
		refresh_M(); //to refresh alive no_alive 
	}

	bool Compart::split(){
		if(reps.size() >= (unsigned int) par_splitfrom) {
			Compart *target = parent->get(); //choose a random cell

			if(target == this){ // in case it should replicate to its own position it only looses half of its content
				//for(auto &rep : reps){
				for(std::list<Compart::ScmRep*>::iterator rep = reps.begin(); rep != reps.end(); ){
					Compart::ScmRep *temp = *(rep++);
					if(gsl_rng_uniform(r) < 0.5) die(temp);
				}
			} else { // overwriting other vesicules
				target->clear(); //kill cell 
				
				// hypergeometric
				for(std::list<Compart::ScmRep*>::iterator emigrant = reps.begin(), temp_it; emigrant != reps.end(); ){
					temp_it = emigrant++; //it is neccessary to keep emigant in the range of reps
					if(gsl_rng_uniform(r) < 0.5) {
						(**temp_it).assignCompart(target);
					} 
				}
				target->refresh_M(); // to refresh no_alive
			}
			refresh_M(); // to refresh no_alive
			
			++(parent->no_last_splits);
			return true;
		}

		return false;

	}

	void Compart::ScmRep::replicates(){
		vesicule->replicate(this);
	}

	void Compart::ScmRep::degradets(){
		if(!empty) die(); // make it an empty replicator
		auto oldves = assignCompart(nullptr);
		//updateM(); // refresh M() and refresh everyones Rep (no need to refresh brokenStick here)
		oldves->refresh_M(); // refresh M() and refresh everyones Rep (no need to refresh brokenStick here)
		updateDeg(); // let brokenStick know about the loss
	}

	void Compart::die(Compart::ScmRep *rep){
		if(!rep->empty) rep->die(); // make it an empty replicator
		rep->updateDeg(); // let brokenStick know about the loss
		rep->updateRep(0.0); // let brokenStick know about the loss
		rep->assignCompart(nullptr);
	}

//	std::list<Compart::ScmRep>::iterator Compart::add(std::list<Compart::ScmRep>::iterator it, std::list<Compart::ScmRep> &from){
//		reps.splice(reps.begin(), from, it);
//		return(reps.begin());
//	}

	std::list<Compart::ScmRep*>::iterator Compart::add(Compart::ScmRep *rep){
		reps.push_front(rep);
		return(reps.begin());
	}

	Compart::ScmRep* Compart::add(){
		//reps.push_front(parent->rep_stack.back()); // get empty replicator
		Compart::ScmRep* newrep = parent->rep_stack.back(); // get empty replicator
		parent->rep_stack.pop_back(); // delete from stack
					      
		newrep->assignCompart(this); // let it know where it has gotten to		

		return( newrep );
	}
	
	Compart::ScmRep* Compart::add(std::string newseq){
		
		auto newrep = add();

		//add new value
		//*dynamic_cast<rnarep::CellContent*>(&*newrep) = newseq;
		*static_cast<rnarep::CellContent*>(newrep) = newseq;
		

		if ( newrep->seq.length() == 0) {
			std::cerr << "Warning: Compart::add(string): invalid sequence " << newseq << " !" << std::endl;
			die( newrep );
		}

		return( newrep );
	}

	inline bool Compart::alive() {
		return _alive;
	}

	// M() is suposed to be called every time when changed
	void Compart::refresh_M(){
		// calculate M
		M = 1;

		if(reps.empty()){
			M = 0;
		} else {
			for(int a = 0; a < par_noEA; a++){
				double akt = 0;
				for(auto &rep : reps){
					// M(x) = prod(sum (a_i))
					akt += rep->geta(a) ;
				}
				M *= akt;
			}

			M = std::pow(M, reciproc_noEA);
		}
		
		//update alive and no_alive
		if(_alive != (bool) M){ //it has changed and...
			if((bool) M){ // ...new state is alive
				no_alive++;
				_alive=true;
			} else { // ...new state is dead
				no_alive--;
				_alive=false;
			}
		}

		for(auto &rep : reps) rep->updateRep(M);
	}

	inline double Compart::get_M() const {return M;} 


	void Compart::replicate(Compart::ScmRep* const templ){
			//replicate it!
			++(parent->no_last_replicates);
			auto newrep = add(); //now reps.begin() = newrep
			newrep->replicate( *static_cast<rnarep::CellContent*>(templ) );

			newrep->updateDeg(); // refresh degr broken stick of child - repl brokenstick will be updated anyway

			//splitting
			if(!split()){
				refresh_M(); // in case of splitting split takes care of refreshing
			}
	}

	/*void Compart::update(){
		//eplication, degradation and splitting only if cell is not empty
		if(reps.size()){
			//replication
			auto degr_from = replication();

			//DEGRADATION
			for(std::list<Compart::ScmRep*>::iterator rep = degr_from; rep != reps.end(); ){
				Compart::ScmRep *temp_rep = *(rep++); //to keep rep in the range of reps 
				if(temp_rep->Pdeg > gsl_rng_uniform(r) ) {
					++(parent->no_last_deaths);
					die(temp_rep);
				}
			}

			//splitting
			if( !split() ) get_M(); // in case of splitting no_alive is refreshed, if it did not happen, we should refresh it, as DEG may have changed it
		}
	}*/

	unsigned int CompartPool::discoverComparts(const char * sourcedir){
		std::string command;
		struct stat sb;

		// test dir
		if(stat(sourcedir, &sb) != 0 || !(sb.st_mode & S_IFDIR)) return 0; // if sourcedir does not exist or is a file it returns 0

		//iterate tru file names
		unsigned int no_files = 0;
		for (const auto & file : fs::directory_iterator(sourcedir)){
			if(file.is_regular_file()){
				std::stringstream fn(file.path().stem().string());
				std::string token;

				// read in and check first word
				getline(fn, token, '_');
				if(token != "bubble") continue;

				// read in second word (time) and check first letter
				getline(fn, token, '_');
				if(token[0] != 't') continue;

				// read in time
			 	unsigned int at = std::stoi(token.substr(1));
				if( bubblefiles.find(at) == bubblefiles.end() ) { // if no file has been given at this timepoint before...
					bubblefiles.insert( Bubbles::value_type( at , std::list<std::string>() ) ); // ... then init it will null value
				}

				bubblefiles.find(at)->second.push_back( file.path() );

				//std::cout << std::stoi(token.substr(1))  << std::endl;
				++no_files;
			}
		}

		return no_files;
	}

	void CompartPool::autoCompartInput(){
		auto match = bubblefiles.find(time);
		if(match != bubblefiles.end()){ // if found
			unsigned int cellno = 0;
			for(auto input = match->second.begin(), end = match->second.end(); input != end; ++input){
				compartFromFile(input->c_str(), cellno++);
			}
		}
	}

	bool CompartPool::compartFromFile(const char * infile){
		std::string line, word;
		std::ifstream file(infile);

		if(!file.is_open()) {
			std::cerr << "ERROR: compartFromFile: file can not be opened!" << std::endl;
			return false;
		}

		//get a random cell and empty it
		auto target = get();
		target->clear();

		while(std::getline(file, line)){
			std::istringstream linestream(line);
			linestream >> word;
			(target->add(word))->updateDeg();
		}

		target->refresh_M();

		return true;
	}

	bool CompartPool::compartFromFile(const char * infile, const unsigned int n){
		if(n>=size) throw std::invalid_argument("n exceeds size!");

		std::string line, word;
		std::ifstream file(infile);

		if(!file.is_open()) {
			std::cerr << "ERROR: compartFromFile: file can not be opened!" << std::endl;
			return false;
		}

		//get a cell and empty it
		Compart *target = get(n);
		target->clear();

		while(std::getline(file, line)){
			std::istringstream linestream(line);
			linestream >> word;
			(target->add(word))->updateDeg();
		}

		target->refresh_M();

		return true;
	}

	// Constructor
	CompartPool::CompartPool(int _size):
		degpool(size*(par_splitfrom-1)+2),
		reppool(size*(par_splitfrom-1)+2),
		size(_size),
		no_last_splits(0),
		no_last_replicates(0),
		no_last_deaths(0),
		no_reps_last(0),
		time(0)
	{
		comparts = new class Compart* [size];
		for(Compart **comp = comparts, **endcomp = comparts+size; comp != endcomp; comp++){
			*comp = new class Compart;
			(*comp)->parent = this;
		}

		replicators = new class Compart::ScmRep[size*(par_splitfrom-1)+1];

		degpool.push_back(1,nullptr); // first element is no_replicators-sum(Pdeg)
		reppool.push_back(1, nullptr); // first element is no_replicators * Cnorep

		for(Compart::ScmRep *ptr = replicators, *end = replicators + (size*(par_splitfrom-1) +1); ptr != end; ++ptr) {
			rep_stack.push_back(ptr);
			degpool.push_back(0,ptr);
			reppool.push_back(0,ptr);
			ptr->setBindings(&reppool.back(), &degpool.back());
		}

//		std::cout << "Basic Constructor Called" << std::endl;
	}

	void CompartPool::init_fromfile(char *infile) {
		std::string line, word;

		std::ifstream file(infile);

		if(!file.is_open()) std::cerr << "ERROR: init_fromfile: file can not be opened!" << std::endl;

		Compart::ScmRep::no_replicators = 0; //clearing number of replicators 

		for(unsigned int cnum = 0; cnum < size; ++cnum){
			for(int no_input = par_num_input_content; no_input-- && std::getline(file, line); ){
				std::istringstream linestream(line);
				linestream >> word;
				comparts[cnum]->add(word);
			}
		}

	}

	CompartPool::~CompartPool(){
		if(size) {
			for(unsigned int i = 0; i < size; ++i ){
				delete comparts[i];
			}
			delete [] (comparts);
			delete [] (replicators);
			//delete [] (temp_comparts);
		}
//				std::cout << "Deconstructor Called" << std::endl;
	}

	inline Compart* CompartPool::get(int cell) {
		return( *(comparts + cell) );
	}

	///gives back pointer to random comp
	inline Compart* CompartPool::get() {
		return( *(comparts + gsl_rng_uniform_int(r, size)) );
	}

	inline Compart* CompartPool::get(Compart *except) {
		if(size < 2) return NULL;

		int n = (int) (except - *comparts);
		return( *(comparts + (n + gsl_rng_uniform_int(r, size-1) + 1) % size) );

		/*compart* out = comparts + gsl_rng_uniform_int(r, size);
		while(except == out) comparts + gsl_rng_uniform_int(r, size);
		return(out);*/
	}

	///Random update
	/*int CompartPool::rUpdate(int gens){
		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		int mtime = 0;
		for(mtime = time + gens; time < mtime; time++){ //updating generations
			autoCompartInput();

			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			//UPDATING
			no_last_splits = no_last_replicates = no_last_deaths = 0;
			no_reps_last = Compart::ScmRep::no_replicators;
			for(unsigned int iter = size+1; --iter; ) comparts[ gsl_rng_uniform_int(r, size) ]->update();

			// quit conditions
			if(testQuit()) break;

		}
		
		//saving/outputting
		if(par_output_interval) do_output();
		if(par_save_interval) if(save()) return -2;


		if(time < mtime) return(2); // quit condition triggered
		if(Compart::ScmRep::no_replicators) return(1); // it has died out
		else return(0);
	}*/

	/// Custom update for this simulation
	int CompartPool::cUpdate(int gens){
		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		int mtime = 0;
		for(mtime = time + gens; time < mtime; time++){ //updating generations
			autoCompartInput();

			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();
			no_last_splits = no_last_replicates = no_last_deaths = 0; // has to be after writing output
			no_reps_last = Compart::ScmRep::no_replicators;

			//UPDATING
			for(unsigned int iter = Compart::ScmRep::no_replicators; --iter; ) {
				// replicate
				reppool[0].update(rnarep::CellContent::no_replicators * par_claimNorep); // the claim of not happening replication is no_replicators * Claim_norep
				double rn = gsl_rng_uniform(r);
				auto target = reppool.draw(rn);
				if(target != nullptr) target->replicates();

				// degrade
				degpool[0].update(rnarep::CellContent::no_replicators -  degpool.cumsum() + degpool[0].get_p() ); // the claim of not happening degradation is no_replicators - sum(Pdeg[1-N])
				target = degpool.draw(gsl_rng_uniform(r));
				if(target != nullptr) target->degradets();
			}

			// quit conditions
			if(testQuit()) break;

		}
		
		//saving/outputting
		if(par_output_interval) do_output();
		if(par_save_interval) if(save()) return -2;


		if(time < mtime) return(2); // quit condition triggered
		if(Compart::ScmRep::no_replicators) return(1); // it has died out
		else return(0);
	}

	//Update according to a random order (in every generation all cells will be updated)
	/*int CompartPool::Update(int gens){
		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		int mtime=0; // so it can detect quit condition 
		for(mtime = time + gens ; time < mtime ; time++){ //updating generations
			autoCompartInput();

			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			//UPDATING
			for(unsigned int iter = 0; iter < size; iter++) comparts[ iter ]->update();

			// quit conditions
			if(testQuit()) break;

		}

		if(par_save_interval) save();
		if( std::strlen(par_output_filename) ) do_output();


		if(time < mtime) return(2); // quit condition triggered
		if(Compart::ScmRep::no_replicators) return(1); // it has died out
		else return(0);
	}*/

	//Update according to a random order (in every generation all cells will be updated)
	/*int CompartPool::oUpdate(int gens){
		int *order;
		unsigned int iter=0, temp = 0, target = 0;

		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}


		//init order
		order = new int[size];
		for(iter = 0; iter < size; iter++){
			order[iter] = iter;
		}

		//for(int mtime = time + gens ; time < mtime && (Compart::no_alive || Compart::ScmRep::no_replicators) ; time++){ //updating generations
		int mtime=0; // so it can detect quit condition 
		for(mtime = time + gens ; time < mtime ; time++){ //updating generations
			autoCompartInput();

			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			//UPDATING
			for(iter = 0; iter < size; iter++){
				target = gsl_rng_uniform_int(r, size - iter);
				if (target) { //it changes with a value higher that itself
					target += iter; //correct it to point to the absolute position
					//switch target and iter
					temp = order[target];
					order[ target ] = order[ iter ];
					order[ iter ] = temp;
					//updateStep( temp );
					comparts[ temp ]->update(); //temp is equaal to order[iter]
				}
				else { //it would change with itself so no change at all
					//updateStep( order[iter] );
					comparts[ order[iter] ]->update();
				}
			}

			// quit conditions
			if(testQuit()) break;

		}

		if(par_save_interval) save();
		if( std::strlen(par_output_filename) ) do_output();

		delete [] (order);

		if(time < mtime) return(2); // quit condition triggered
		if(Compart::ScmRep::no_replicators) return(1); // it has died out
		else return(0);
	}*/
	
	bool CompartPool::testQuit() const {
			if(!par_quit) return false;

			if(par_quit & qreplicator) if(Compart::ScmRep::no_replicators == 0) return true;
			if(par_quit & qcompart) if(Compart::no_alive == 0) return true;
			if(par_quit & qsplit) if(no_last_splits > 0) return true; 
			if(par_quit & qfull) {
				unsigned int comp = 0;
				while( comp < size && !(comparts[comp]->reps.empty()) ) {++comp;} 
				if(comp == size) return true;
			}
			if(par_quit & qalive) if(Compart::no_alive == size) return true;

			return false;

	}

	int CompartPool::openOutputs(){
		std::string name, command;
		 
		name += par_outdir;
		
		//create directory output if it does not exist (Linux only!!)
		command += "mkdir -p ";
		command += par_outdir;
//		std::cout << command << std::endl;
		system(command.c_str());

		//create directory for output
		name += "/";
		name += par_ID;

		command.clear();
		command += "test -d ";
		command += name;
		
//		std::cout << command << std::endl;

		if(system(command.c_str())) { //directory does not exist
//			std::cout << "not exists" << std::endl;

			command.clear();
			command += "mkdir ";
			command += name;
			
//			std::cout << command << std::endl;
			system(command.c_str());
		}
		else{ //directory already exist
//			std::cerr << "exists" << std::endl;
			//check if files already exists
			command.clear();
			command += "test -f ";
			command += name;
			command += "/";
			command += par_output_filename;
//			std::cout << command << std::endl;
			if(!system(command.c_str())){ //exists already
				//try a new directory name
				name += "_";
				name += std::to_string(gsl_rng_uniform_int(r, 1000));
				command.clear();
				command += "test -d ";
				command += name;
//				std::cout << command << std::endl;
				if(!system(command.c_str())) { //already exzist
					std::cerr << "ERROR: conflicting simulations with identical IDs. Quitting..." << std::endl;
					return (1);
				}
				else{ // new dir does not exist
					//std::cerr << "WARNING: directory already existed with output files! Tried a new name: " << name << std::endl;
					command.clear();
					command += "mkdir ";
					command += name;
					
//					std::cout << command << std::endl;
					system(command.c_str());
				}
			} //files exist
		} //directory exists

		name += "/";
		//output directory exists and its path is in name

		//creating output file
		command.clear();
		command += name;
		command += par_output_filename;

		output.open(command);
		if(!output.is_open()) {
			std::cerr << "ERROR: output file (" << name << par_output_filename << ") cant be opened" << std::endl;
			return (2);
		}

		//adding header to output
		output << "time;replicators;no_alive;no_last_splits;mean_M;sd_M";
		output << ";no_par;mean_R_par;mean_length_par;mean_mfe_par" ;
		for(int e = 0; e < par_noEA; e++){
			output << ";no_enz" << e << ";mean_R_enz" << e << ";mean_length_enz" << e << ";mean_mfe_enz" << e << ";mean_a_enz" << e ;
		}

		for(int a = 0; a <= par_noEA; a++) {
			output << ";no_A" << a ;
		}

		output << ";percent_replicated;percent_died" << std::endl;

		output.flush();

		//prepare vectors for output
		
		out_no.reserve(par_noEA + 1);
		out_noA.reserve(par_noEA + 1);
		out_R.reserve(par_noEA + 1);
		out_length.reserve(par_noEA + 1);
		out_a.reserve(par_noEA + 1); //parasites do not have activities, no need to calculate with 0th value!
		out_mfe.reserve(par_noEA + 1);

		//creating SAVE directory
		command.clear();
		command += "mkdir -p ";
		
		name += par_savedir;

		command += name;
//		std::cout << command << std::endl;
		system(command.c_str());

		savedir = name;

		return 0;

	}

	void CompartPool::do_output(){

//		std::cout << "output" << std::endl;
		/* what i need:
			time,								#ok
			no_replicators,							#ok
			alive,								#added
			mean_M,								#added
			sd_M,								#added
			[by akt: No, Rs mean, length mean, mfe mean, alpha mean], 	#ok
			[by no akts: number]						#ok

		*/

		//clearing
		out_no.assign(par_noEA + 1, 0);
		out_noA.assign(par_noEA + 1, 0);
		out_R.assign(par_noEA + 1, 0);
		out_length.assign(par_noEA + 1, 0);
		out_a.assign(par_noEA + 1, 0);
		out_mfe.assign(par_noEA + 1, 0);

		//compart level variables
		double sum_M = 0;
		double sum_M2 = 0;

		//calculating values
		for(Compart **acell = comparts, **end = (Compart **) comparts + size ; acell != end; acell++){
			Compart *cell = *acell;
//			std::cout << "examined cell" << std::endl;
			
			//compart level calculations
			if(cell->alive()){
				double metab = cell->get_M();
				sum_M += metab;
				sum_M2 += metab*metab;
			}

			//replicators in compart
			for(auto &rep : cell->reps){ // iterating tru reps in cell
				//how much activities does it have?
				int no_acts = rep->get_no_acts();
//				if(no_acts >= out_noA.size()){
//					std::cerr << "Warning: do_output: something went wrong!" << std::endl;
//				}
				out_noA[no_acts]++;

				if(no_acts){ //it is not a parasite
					for(int ea = 1; ea <= par_noEA; ea++){
						double activity = rep->geta(ea - 1); //indexing of "out_" arrays starts at parazite, "a" starts with activity 0
//						std::cout << "in rep (" << cell->vals->get_type() <<  ")  having " << no_acts << " activities " << ea << "th activity is " << activity << std::endl;
						if(activity) { //if it has activity ea
							out_no[ea]++;
							out_R[ea] += rep->getR();
							out_length[ea] += rep->get_length();
							out_mfe[ea] += rep->get_mfe();
							out_a[ea] += activity;
						}
					}
				}
				else { //it is a parasite
					out_no[0]++;
					out_R[0] += rep->getR();
					out_length[0] += rep->get_length();
					out_mfe[0] += rep->get_mfe();
				} 

				//calculating means
				/*for(int ea = 0; ea <= par_noEA; ea++){
					out_R[ea] /= out_no[ea];
					out_length[ea] /= out_no[ea];
					out_a[ea] /= out_no[ea];
					out_mfe[ea] /= out_no[ea];
				}*/

			} // for reps
		} // tru cells in matrix

		//outputting
		///CompartPool and Compart level variables
		output << time << ';' 						//time
			<< Compart::ScmRep::no_replicators << ';' 		//no_replicators
			<< Compart::no_alive << ';' 				//no_alive
			<< no_last_splits << ';' 				
			<< (Compart::no_alive?(sum_M/(double)Compart::no_alive):0) << ';' 		//mean_M
			<< dvtools::sd(Compart::no_alive, sum_M, sum_M2);	//sd_M
		///Replicator level variables
		double no;
		for(int ea = 0; ea <= par_noEA; ea++) {
			if( (no = (double) out_no[ea]) ){
//				std::cout << "outputting: no= " << no << std::endl;
				output << ';' << no << ';' << out_R[ea]/no << ';' << out_length[ea]/no << ';' << out_mfe[ea]/no;
				if(ea) output << ';' << out_a[ea]/no;
			}
			else {
//				std::cout << "outputting: no= " << no << " (supposedly 0)" << std::endl;
				output << ";0;0;0;0";
				if(ea) output << ";0";
			}

		}

		for(int ea = 0; ea <= par_noEA; ea++) {
			output << ';' << out_noA[ea] ;
		}

		output << ';' << static_cast<double>(no_last_replicates) << ';' << static_cast<double>(no_last_deaths)/no_reps_last << std::endl;

		output.flush();

	}


	int CompartPool::save(){
		/* Outputs:
		- XML file containing EVERYTHING (except parameters)
		- rng binary state file
		*/

		std::string filename;

		//1. SAVING CONTENT
		{
			//open outputs
			if(!savedir.length()) {
				std::cerr << "ERROR: No savedir inicialised! Please do run CompartPool::openOutputs() before saving!" << std::endl;
				return (1);
			}
			
			filename = savedir; 
			filename += '/'; 
			filename += std::to_string(time);
			filename += ".xml";
			std::ofstream out(filename);
			
			if(!out){
				std::cerr << "ERROR: Could not open file for saving content!" << std::endl;
				return 2;
			}

			//saving content
			assert(out.good());
			boost::archive::xml_oarchive oa(out);
			oa << boost::serialization::make_nvp("mcrscm", *this);

		}

		//2. saving rng state
		filename.clear();
		filename = savedir;
		filename += "/rngstate";
		filename += std::to_string(time);
		filename += ".bin";
		//std::ofstream rngout(filename, std::ios::out | std::ios::binary);

		if(randomszam_mentes(filename.c_str(), r)) {
			std::cerr << "ERROR: Could not open file for saving random number generator state!" << std::endl;
			return 3;
		}

		return 0;
	}

} // namespace cmdv


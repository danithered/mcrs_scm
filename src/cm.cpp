#include "cm.h"
#include "include/rnarep.h"

namespace cmdv {
	//int no_births=0;
	//int no_deaths=0;
	
	unsigned int Compart::no_alive = 0;

	std::list<rnarep::CellContent> Compart::wastebin;

	Compart::Compart(): parent(NULL), reciproc_noEA(1 / (double) par_noEA), _M(0.0), _alive(false), changed_content(true){}

	/*void Compart::operator =(Compart& origin){
		clear();
		//for(auto rep = origin.reps.begin(); rep != origin.reps.end(); rep++){
			reps.assign( origin.reps.begin(), origin.reps.end() ); 
		//}
		//no need to delete temp_compart, split() will do that anyway...
	}*/

	void Compart::clear(){
		for(auto rep = reps.begin(); rep != reps.end(); ) die(rep++);
		M(); //to refresh alive no_alive and no_replicators according to awake
	}

	bool Compart::split(){
		if(reps.size() > (unsigned int) par_splitfrom) {
			Compart *target = parent->get(); //choose a random cell

			if(target == this){ // in case it should replicate to its own position it only looses half of its content
				//for(auto &rep : reps){
				for(std::list<rnarep::CellContent>::iterator rep = reps.begin(), temp_it; rep != reps.end(); ){
					temp_it = rep++;
					if(gsl_rng_uniform(r) < 0.5) die(temp_it);
				}
			} else { // overwriting other vesicules
				target->clear(); //kill cell 
				auto &targetreps = target->reps;
				
				// hypergeometric
				for(std::list<rnarep::CellContent>::iterator emigrant = reps.begin(), temp_it; emigrant != reps.end(); ){
					temp_it = emigrant++; //it is neccessary to keep emigant in the range of reps
					if(gsl_rng_uniform(r) < 0.5) {
						targetreps.splice(targetreps.begin(), reps, temp_it); //splice puts an elment of a list to an other. Args: iterator target pos, origin list, origin iterator to be transferred
					} 
				}
				changed_content = true; // target is flagged with changed content in clear()-s die()
				target->M(); // to refresh no_alive
			}
			M(); // to refresh no_alive
			
			++(parent->no_last_splits);
			return true;
		}

		return false;

	}

	void Compart::die(std::list<rnarep::CellContent>::iterator rep){
		if(!rep->empty) rep->die();
		wastebin.splice(wastebin.begin(), reps, rep);
		changed_content = true;
	}

	std::list<rnarep::CellContent>::iterator Compart::add(std::list<rnarep::CellContent>::iterator it, std::list<rnarep::CellContent> &from){
		reps.splice(reps.begin(), from, it);
		changed_content=true;
		return(reps.begin());
	}

	std::list<rnarep::CellContent>::iterator Compart::add(rnarep::CellContent rep){
		reps.push_front(rep);
		changed_content=true;
		return(reps.begin());
	}

	std::list<rnarep::CellContent>::iterator Compart::add(){
		if(wastebin.empty()){ //have to create new
			reps.emplace_front();
		} else { //reuse some from wastebin
			reps.splice(reps.begin(), wastebin, wastebin.begin() );
		}
		changed_content=true;

		return( reps.begin() );
		
	}
	
	std::list<rnarep::CellContent>::iterator Compart::add(std::string newseq){
		
		auto newrep = add();

		//add new value
		*newrep = newseq;
		
		if ( newrep->seq.length() == 0) {
			std::cerr << "Warning: Compart::add(string): invalid sequence " << newseq << " !" << std::endl;
			die( newrep );
		}

		return( newrep );
	}

	inline bool Compart::alive() {
		if (changed_content) M();
		return _alive;
	}

	double Compart::M(){
		if(changed_content){
			// calculate M
			_M = 1;

			if(reps.empty()){
				_M = 0;
			} else {
				for(int a = 0; a < par_noEA; a++){
					double akt = 0;
					for(auto met = reps.begin(); met != reps.end(); met++){
						// M(x) = prod(sum (a_i))
						akt += met->geta(a) ;
					}
					_M *= akt;
				}

				_M = std::pow(_M, reciproc_noEA);
			}
			
			//update alive and no_alive
			if(_alive != (bool) _M){ //it has changed and...
				if((bool) _M){ // ...new state is alive
					no_alive++;
					_alive=true;
				} else { // ...new state is dead
					no_alive--;
					_alive=false;
				}
			}

			// flag it as calculated
			changed_content = false;

		}

		return _M;
	}

	std::list<rnarep::CellContent>::iterator Compart::replication(){
		auto oldfirst= reps.begin();
		
		if (const double met = M(); met != 0.0 ) { // is alive
			(parent->no_reps_last_in_alive) += reps.size();
			const double CperM = par_claimNorep / met;
			for(auto &rep : reps){
				const double replicability = rep.getR();
				if(gsl_rng_uniform(r) < (replicability / (replicability + CperM)) ){
					//replicate it!
					++(parent->no_last_replicates);
					auto newrep = add(); //now reps.begin() = newrep
					newrep->replicate( rep );
				}
			}
		}

		return(oldfirst);
	}

	/// replication method with predetermined metabolism. Return value is if it is above splitsize.
	bool Compart::replication(const double met){
		if(met == 0.0) return false; // if there is no metabolism there can be no replication, therefor (if everything went ok) can be no splitting

		// make order map
		std::map<const double, rnarep::CellContent *> reps_by_R;

		for ()


			(parent->no_reps_last_in_alive) += reps.size();
			const double CperM = par_claimNorep / met;
			for(auto &rep : reps){
				const double replicability = rep.getR();
				if(gsl_rng_uniform(r) < (replicability / (replicability + CperM)) ){
					//replicate it!
					++(parent->no_last_replicates);
					auto newrep = add(); //now reps.begin() = newrep
					newrep->replicate( rep );
				}
			}

		return(oldfirst);
	}

	void Compart::update(){
		//replication, degradation and splitting only if cell is not empty
		if(reps.size()){
			// Metabolism
			const double metab = M();

			//DEGRADATION
			while(std::list<rnarep::CellContent>::iterator rep = reps.begin(), temp_rep; rep != reps.end()){
				temp_rep = rep++; //to keep rep in the range of reps 
				if(temp_rep->Pdeg > gsl_rng_uniform(r) ) {
					++(parent->no_last_deaths);
					die(temp_rep);
				}
			}

			//replication
			replication_withM(metab);

				//splitting
				if( !split() ) M(); // in case of splitting no_alive is refreshed, if it did not happen, we should refresh it, as DEG may have changed it
		}
	}

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
			target->add(word);
		}


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
			target->add(word);
		}


		return true;
	}

	// Constructor
	CompartPool::CompartPool(int _size): size(_size), no_last_splits(0), no_last_replicates(0), no_last_deaths(0), no_reps_last_in_alive(0), no_reps_last(0){
		time=0;
		//saving_freq = 0;

		comparts = new class Compart* [size];
		for(Compart **comp = comparts, **endcomp = comparts+size; comp != endcomp; comp++){
			*comp = new class Compart;
			(*comp)->parent = this;
		}

		//no_replicators = 0;

//				std::cout << "Basic Constructor Called" << std::endl;
	}

	void CompartPool::init_fromfile(char *infile) {
		std::string line, word;

		std::ifstream file(infile);

		if(!file.is_open()) std::cerr << "ERROR: init_fromfile: file can not be opened!" << std::endl;

		rnarep::CellContent::no_replicators = 0; //clearing number of replicators 

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
	int CompartPool::rUpdate(int gens){
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
			no_reps_last_in_alive = no_last_splits = no_last_replicates = no_last_deaths = 0;
			no_reps_last = rnarep::CellContent::no_replicators;
			for(unsigned int iter = size+1; --iter; ) comparts[ gsl_rng_uniform_int(r, size) ]->update();

			// quit conditions
			if(testQuit()) break;

		}
		
		//saving/outputting
		if(par_output_interval) do_output();
		if(par_save_interval) if(save()) return -2;


		if(time < mtime) return(2); // quit condition triggered
		if(rnarep::CellContent::no_replicators) return(1); // it has died out
		else return(0);
	}

	//Update according to a random order (in every generation all cells will be updated)
	int CompartPool::Update(int gens){
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
		if(rnarep::CellContent::no_replicators) return(1); // it has died out
		else return(0);
	}

	//Update according to a random order (in every generation all cells will be updated)
	int CompartPool::oUpdate(int gens){
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

		//for(int mtime = time + gens ; time < mtime && (Compart::no_alive || rnarep::CellContent::no_replicators) ; time++){ //updating generations
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
		if(rnarep::CellContent::no_replicators) return(1); // it has died out
		else return(0);
	}
	
	bool CompartPool::testQuit() const {
			if(!par_quit) return false;

			if(par_quit & qreplicator) if(rnarep::CellContent::no_replicators == 0) return true;
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
				double metab = cell->M();
				sum_M += metab;
				sum_M2 += metab*metab;
			}

			//replicators in compart
			for(auto rep = cell->reps.begin(); rep != cell->reps.end(); rep++){ // iterating tru reps in cell
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
			<< rnarep::CellContent::no_replicators << ';' 		//no_replicators
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

		output << ';' << static_cast<double>(no_last_replicates)/no_reps_last_in_alive << ';' << static_cast<double>(no_last_deaths)/no_reps_last << std::endl;

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


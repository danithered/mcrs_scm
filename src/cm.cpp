#include "cm.h"

namespace cmdv {
	//int no_births=0;
	//int no_deaths=0;
	
	unsigned int Compart::no_alive;

	void Compart::clear(){
		for(auto rep = reps.begin(); rep != reps.end(); rep++) die(rep);
		leftover = 0; //so it did not inherits lost compartments metabolism
	}

	int Compart::split(){
		//Compart *target = parent->get(this); //choose a random cell
		//give it to temp
		parent->temp_comparts.emplace_back();
		Compart *target = &(parent->temp_comparts.back());
						     
		std::list<rnarep::CellContent> *targetreps = &(target->reps);

		target->clear(); //kill cell
		
		//hypergeometric
		for(auto emigrant = reps.begin(); emigrant != reps.end(); ){
			if(gsl_rng_uniform(r) < 0.5) {
				targetreps->splice(targetreps->begin(), reps, emigrant++);
			} 
			else emigrant++;
		}
		
		//target->updateable = false;
		
		/*clear leftover
		 * it would be unfair to inherit something from a metabolism with a possibly completely different composition 
		 * consider it a price for splitting */
		leftover = 0;
		//target->leftover = 0; it is handled in clear()

		return(0);

	}

	void Compart::die(std::list<rnarep::CellContent>::iterator rep){
		if(!rep->empty) rep->die();
		wastebin.splice(wastebin.begin(), reps, rep);
		
	}

	std::list<rnarep::CellContent>::iterator Compart::add(std::list<rnarep::CellContent>::iterator it, std::list<rnarep::CellContent> &from){
		reps.splice(reps.begin(), from, it);
		return(reps.begin());
	}

	std::list<rnarep::CellContent>::iterator Compart::add(rnarep::CellContent rep){
		reps.push_front(rep);
		return(reps.begin());
	}

	std::list<rnarep::CellContent>::iterator Compart::add(){
		if(wastebin.empty()){ //have to create new
			reps.emplace_front();
		} else { //reuse some from wastebin
			reps.splice(reps.begin(), wastebin, wastebin.begin() );
		}

		return( reps.begin() );
		
	}
	
	std::list<rnarep::CellContent>::iterator Compart::add(std::string newseq){
		
		if(wastebin.empty()){ //have to create new
			reps.emplace_front();
		} else { //reuse some from wastebin
			reps.splice(reps.begin(), wastebin, wastebin.begin() );
		}

		//add new value
		reps.front() = newseq;

		//annotate
		//reps.front().annotate(); // operator= handles it itself 

		return(reps.begin());
	}

	double Compart::M(){
		double M = 1;

		for(int a = 0; a < par_noEA; a++){
			double akt = 0;
			for(auto met = reps.begin(); met != reps.end(); met++){
				// M(x) = prod(sum (a_i))
				akt += met->geta(a) ;
			}
			M *= akt;
		}
		
		//update alive and no_alive
		if(alive != (bool) M){ //it has changed
			if(M){ //new state is alive
				alive=true;
				no_alive++;
			} else { //new state is dead
				alive=false;
				no_alive--;
			}
		}

		return std::pow(M, reciproc_noEA);
	}

	std::list<rnarep::CellContent>::iterator Compart::replication(){
		auto oldfirst= reps.begin();
		
		//compute number of new replicators
		leftover += M() * par_MN;  //number of new replicators is linear function of M, slope is par_MN, intercept is 0 (+ leftover from previous round)
		int number_of_new = leftover; //number is an integer, float part will be used in next round
		leftover -= number_of_new; //leftover contains that part of metabolism that is not used in this round

		if(number_of_new){
			double sum=0;
			double *replicabilities = new double [reps.size()];
			double *r_it = replicabilities;
			int size_origin = reps.size();

			//extract Rs
			for(auto rep = reps.begin(); rep != reps.end(); rep++) { 
				if(!rep->empty){
					sum += (*r_it = rep->getR()) ;
				}
				r_it++;
			}

			//create new replicators
			for(int already_added = 0; already_added < number_of_new; already_added++){
				//which replicator is being replicated
				int decision = dvtools::brokenStickVals(replicabilities, size_origin, sum, gsl_rng_uniform(r)) ; //It can be negative! See: brokenStickVals
				//replicate it!
				auto newrep = add();
				//newrep->replicate( reps[decision + already_added] );
				newrep->replicate( *std::next(reps.begin(), decision + already_added) );
			}
			

			//free
			delete [] (replicabilities);
		}

		return(oldfirst);

	}

	void Compart::update(){
		if(updateable){
			//replication and degradation only if cell is not empty
			if(reps.size()){
				//replication
				auto degr_from = replication();
				
				//DEGRADATION
				for(auto rep = degr_from; rep != reps.end(); rep++){
					if(rep->Pdeg > gsl_rng_uniform(r) ) {
//						no_deaths++;
						die(rep);
					}
				}
			}
			//splitting
			if(reps.size() > (unsigned int) par_splitfrom) split();
		}
	}

	///initialises matrix with predefined values, randomly
	/*void CompartPool::init(std::string *pool, double* probs, int no_choices) {
		int i = 0;
		double sum = 0.0;
		
		for(i=0; i < no_choices; i++) {
			sum += probs[i];
		}
		for(i=0; i < size; i++) {
			*(matrix[i].vals) = pool[dvtools::brokenStickVals(probs, no_choices, sum, gsl_rng_uniform(r) )];
		}
	}*/	

	void CompartPool::init_fromfile(char *infile) {
		std::string line, word;

		std::ifstream file(infile);

		if(!file.is_open()) std::cerr << "ERROR: init_fromfile: file can not be opened!" << std::endl;

		rnarep::CellContent::no_replicators = 0; //clearing number of replicators 

		for(int cnum = 0; cnum < size; cnum++){
			for(int no_input = par_num_input_content; no_input-- && std::getline(file, line); ){
				std::istringstream linestream(line);
				linestream >> word;
				comparts[cnum].add(word);
			}
		}

//		std::cout << "Grid initialised with " << rnarep::CellContent::no_replicators << " replicators on a grid of " << size << " cells." << std::endl;
	}

	//a singel update step <- this is called by rUpdate and oUpdate
	/*int CompartPool::updateStep(int cell){
		matrix[cell].update();
		return(0);
	}*/

	void CompartPool::all_updateable(){
		Compart *comp = comparts;
		for(int countdown = size; countdown--; comp++){
			if(!comp->updateable) comp->updateable = true;
		}
	}

	///Random update
	int CompartPool::rUpdate(int gens){
		int iter=0;

		//check if output is open
		if(par_output_interval && !output) {
			std::cerr << "ERROR: output not open" << std::endl;
			return(-1);
		}

		//output/save in case of not init from start
		if(time){
			if(par_output_interval && (time % par_output_interval)) do_output();
			if(par_save_interval && (time % par_save_interval)) save();
		}

		for(int mtime = time + gens ; time < mtime && rnarep::CellContent::no_replicators; time++){ //updating generations
			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			for(iter = 0; iter < size; iter++){
				//UPDATING
				comparts[ gsl_rng_uniform_int(r, size) ].update();
			}
//			std::cout << "Cycle " << time << ": number of total deaths: " << no_deaths << ", number of total births: " << no_births << std::endl;
		}
		
		//saving/outputting
		if(par_output_interval) do_output();
		if(par_save_interval) if(save()) return -2;


		return rnarep::CellContent::no_replicators;
	}

	//Update according to a random order (in every generation all cells will be updated)
	int CompartPool::oUpdate(int gens){
		int *order;
		int iter=0, temp = 0, target = 0;

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

		for(int mtime = time + gens ; rnarep::CellContent::no_replicators && time < mtime ; time++){ //updating generations
			//outputs
			if (par_output_interval && !(time % par_output_interval)) do_output();
			if (par_save_interval && !(time % par_save_interval)) save();

			//make all compart updateable
			all_updateable();

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
					comparts[ temp ].update(); //temp is equaal to order[iter]
				}
				else { //it would change with itself so no change at all
					//updateStep( order[iter] );
					comparts[ order[iter] ].update();
				}
			}
		}

		if(par_save_interval) save();
		if(par_output_filename) do_output();

		delete [] (order);

		return(0);
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
		output << "time;replicators;no_alive;mean_M;sd_M";
		output << ";no_par;mean_R_par;mean_length_par;mean_mfe_par" ;
		for(int e = 0; e < par_noEA; e++){
			output << ";no_enz" << e << ";mean_R_enz" << e << ";mean_length_enz" << e << ";mean_mfe_enz" << e << ";mean_a_enz" << e ;
		}

		for(int a = 0; a <= par_noEA; a++) {
			output << ";no_A" << a ;
		}

		output << std::endl;

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
		for(Compart *cell = comparts, *end = (Compart *) comparts + size ; cell != end; cell++){
//			std::cout << "examined cell" << std::endl;
			
			//compart level calculations
			if(cell->alive){
				double metab = cell->M();
				sum_M += metab;
				sum_M2 += metab*metab;
			}

			//replicators in compart
			for(auto rep = cell->reps.begin(); rep != cell->reps.end(); rep++){ // iterating tru reps in cell
				//how much activities does it have?
				int no_acts = rep->get_no_acts();

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
			<< sum_M/(double)Compart::no_alive << ';' 		//mean_M
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

		output << std::endl;

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
			oa << BOOST_SERIALIZATION_NVP(*this);

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

}


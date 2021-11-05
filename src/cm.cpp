#include "ca.h"

namespace cmdv {
	//int no_births=0;
	//int no_deaths=0;

	int Compart::die(rnarep::CellContent *rep){
		rep->die();
		wastebin.spline(wastebin.begin(), reps, rep);
	}

	rnarep::CellContent* Compart::add(rnarep::CellContent *rep){
		reps.push_front(rep);
		return(rep);
	}

	rnarep::CellContent* Compart::add(){
		if(wastebin.empty()){ //have to create new
			reps.push_front(new rnarep::CellContent);
		} else { //reuse some from wastebin
			reps.push_front( wastebin.front() );
			wastebin.pop_front(); 
		}

		return( reps.front() );
		
	}
	
	rnarep::CellContent* Compart::add(std::string newseq){
		
		if(wastebin.empty()){ //have to create new
			reps.push_front(new rnarep::CellContent);
		} else { //reuse some from wastebin
			reps.push_front( wastebin.front() );
			wastebin.pop_front(); 
		}

		rnarep::CellContent *value = reps.front();

		//add new value
		(*value) = newseq;

		//annotate
		value->annotate();

		return(value);
	}

	double Compart::M(){
		double M = 1, akt = 0;

		for(int a = 0; a < par_noEA; a++){
			akt = 0;
			for(int met = 0; met < no_met_neigh; met++){
				// M(x) = prod(sum (a_i))
				akt += met_neigh[met]->vals->geta(a) ;
			}
			M *= akt;
		}

		return M;
	}

	void Compart::update(){
		double sum = par_claimEmpty;
		int decision = 0;

		if (vals->empty) { // if focal cell is empty (it can be filled with a copy)
			//REPLICATION
			for(int rep = 1; rep < no_repl_neigh; rep++) { // 0. neighbour is self, but it is empty
				if(!repl_neigh[rep]->vals->empty){
					sum += (claims[rep] = repl_neigh[rep]->vals->getR() * repl_neigh[rep]->M() );
				}
				else claims[rep] = 0.0;
			}
			//decision
			if(sum){
				decision = dvtools::brokenStickVals(claims, no_repl_neigh + 1, sum, gsl_rng_uniform(r)) ;
//				if(decision || (sum!=par_claimEmpty) ) {std::cout << "Replication: decision is: " << decision << " from claims:" << std::endl;
//									for(int op = 0; op < no_repl_neigh+1;op++ ) std::cout << claims[op]/sum << "\t";
//									std::cout << std::endl;}
				if(decision){ //claim 0 is claimEmpty NOTE that the probablity of staying empty is not fixed (e.g. 10%)! In case decision is negative see: brokenStickVals
						vals->replicate( *(repl_neigh[decision]->vals) );
						if(gsl_rng_uniform(r) < 0.5) { //havet to switch them at 50 percent
							switchit( *(repl_neigh[decision]) );
						}
//						no_births++;
//						std::cout << "Replication happend. The two molecules:" << std::endl << *(vals->get_seq()) << std::endl << *(repl_neigh[decision]->vals->get_seq()) << std::endl;  
				}
			}
		}
		else { //if focal cell is occupied (it can die)
			//DEGRADATION
			if(vals->Pdeg > gsl_rng_uniform(r) ) {
//				no_deaths++;
//				std::cout << "Degradation with Pdeg " << vals->Pdeg << std::endl;
				vals->die();
			}
		}
	}






	///gives back coordinates of nth cell
	std::vector<int> CompartPool::getCoord(int n) {
		std::vector<int> coords;
		if (layout == square){
			coords.push_back(n % ncol); //x coord
			coords.push_back(n / ncol); //y coord
		}
		else if (layout == hex) {
			//x=n %/% ncol
			//y = n - x*ncol - x%/%2
			//z = 0-x-y
			coords.push_back( n / ncol); //x cubic coord
			coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
			coords.push_back(0 - coords[0] - coords[1]); //z cubic coord
		}

		return(coords);
	}

	///initialises matrix with predefined values, randomly
	void CompartPool::init(std::string *pool, double* probs, int no_choices) {
		int i = 0;
		double sum = 0.0;
		
		for(i=0; i < no_choices; i++) {
			sum += probs[i];
		}
		for(i=0; i < size; i++) {
			*(matrix[i].vals) = pool[dvtools::brokenStickVals(probs, no_choices, sum, gsl_rng_uniform(r) )];
		}
	}	

	void CompartPool::init_fromfile(char *infile) {
		std::string line, word;

		std::ifstream file(infile);
		Compart *cell = matrix;

		if(!file.is_open()) std::cerr << "ERROR: init_fromfile: file can not be opened!" << std::endl;

		rnarep::CellContent::no_replicators = 0; //clearing number of replicators 

		for (int cellnum = 0; cellnum < size && std::getline(file, line); cellnum++){
		//while (std::getline(file, line) ){
			std::istringstream linestream(line);
			linestream >> word;
//			std::cout << "init_fromfile assigning " << word << std::endl;			
			*(cell->vals) = word;
			cell++;
		}

		if(cell != (Compart *) matrix + size) std::cerr << "WARNING: file length is not equal to gridsize!" << std::endl;
//		std::cout << "Grid initialised with " << rnarep::CellContent::no_replicators << " replicators on a grid of " << size << " cells." << std::endl;
	}

	inline Compart* CompartPool::get(int x, int y) {
		if(layout == square){
			return(matrix + y*ncol + x);
		}
		return(NULL);
	}
	///finds cell in pos [x,y] and gives back its number
	inline int CompartPool::getN(int x, int y) {
		if(layout == square){
			return(y*ncol + x);
		}
		else if (layout == hex){
			return(y + x*ncol + x/2);
		}
		return(-1);
	}

	//a singel update step <- this is called by rUpdate and oUpdate
	int CompartPool::updateStep(int cell){
		matrix[cell].update();
		return(0);
	}

	///Random update
	int CompartPool::rUpdate(int gens){
		int iter=0, diff_until = dvtools::fracpart(diff * size * time);

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
				matrix[ gsl_rng_uniform_int(r, size) ].update();
				//DIFFUSION
				for(diff_until += diff; diff_until >= 1; diff_until--){
					matrix[ gsl_rng_uniform_int(r, size) ].diff();
				}
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
		int iter=0, temp = 0, target = 0, diff_until = dvtools::fracpart(diff * time);

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

			//UPDATING
			for(iter = 0; iter < size; iter++){
				target = gsl_rng_uniform_int(r, size - iter);
				if (target) {
					target += iter;
					temp = order[target];
					order[ target ] = order[ iter ];
					order[ iter ] = temp;
					//updateStep( temp );
					matrix[ temp ].update();
				}
				else {
					//updateStep( order[iter] );
					matrix[ order[iter] ].update();
				}
			}
			//DIFFUSION
			for(diff_until += diff; diff_until >= 1; diff_until--){
				for(iter = 0; iter < size; iter++){
					target = gsl_rng_uniform_int(r, size - iter);
					if (target) {
						target += iter;
						temp = order[target];
						order[ target ] = order[ iter ];
						order[ iter ] = temp;
						matrix[ temp ].diff();
					}
					else {
						matrix[ order[iter] ].diff();
					}
				}
			}
		}

		if(par_save_interval) save();
		if(par_output_filename) do_output();

		delete [] (order);

		return(0);
	}

	//gives back coordinates of nth cell
	std::vector<int> CompartPool::getCoord(int n, int type = 0) {
		std::vector<int> coords;
		if (layout == square){
			coords.push_back(n % ncol); //x coord
			coords.push_back(n / ncol); //y coord
		}
		else if (layout == hex) {
			if(type == 0){ //cube coordinates - axial coordinates are the choosen two from this three coords
				coords.push_back( n / ncol); //x cubic coord
				coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
				coords.push_back(0 - coords[0] - coords[1]); //z cubic coord
			}
			else if (type == 1){ //axial coords 
				coords.push_back( n / ncol); //x cubic coord
				coords.push_back( n - coords[0]*ncol - coords[0]/2 ); //y cubic coord
			}
			else if (type == 2){ //offset coords 
				coords.push_back(n % ncol); //x coord
				coords.push_back(n / ncol); //y coord
			}
		}

		return(coords);
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
		output << "time;replicators";
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

	void CompartAut::do_output(){

//		std::cout << "output" << std::endl;
		/* what i need:
			time, alive, [by akt: No, Rs mean, length mean, alpha mean, mfe mean], [by no akts: number]

		*/
		//clearing
		out_no.assign(par_noEA + 1, 0);
		out_noA.assign(par_noEA + 1, 0);
		out_R.assign(par_noEA + 1, 0);
		out_length.assign(par_noEA + 1, 0);
		out_a.assign(par_noEA + 1, 0);
		out_mfe.assign(par_noEA + 1, 0);

		//calculating values
		for(Compart *cell = matrix, *end = (Compart *) matrix + size ; cell != end; cell++){
//			std::cout << "examined cell" << std::endl;
			if(!cell->vals->empty){ // if cell is not empty
				//how much activities does it have?
				int no_acts = cell->vals->get_no_acts();

				out_noA[no_acts]++;

				if(no_acts){ //it is not a parasite
					for(int ea = 1; ea <= par_noEA; ea++){
						double activity = cell->vals->geta(ea - 1); //indexing of "out_" arrays starts at parazite, "a" starts with activity 0
//						std::cout << "in rep (" << cell->vals->get_type() <<  ")  having " << no_acts << " activities " << ea << "th activity is " << activity << std::endl;
						if(activity) { //if it has activity ea
							out_no[ea]++;
							out_R[ea] += cell->vals->getR();
							out_length[ea] += cell->vals->get_length();
							out_mfe[ea] += cell->vals->get_mfe();
							out_a[ea] += activity;
						}
					}
				}
				else { //it is a parasite
					out_no[0]++;
					out_R[0] += cell->vals->getR();
					out_length[0] += cell->vals->get_length();
					out_mfe[0] += cell->vals->get_mfe();
				} 

				//calculating means
				/*for(int ea = 0; ea <= par_noEA; ea++){
					out_R[ea] /= out_no[ea];
					out_length[ea] /= out_no[ea];
					out_a[ea] /= out_no[ea];
					out_mfe[ea] /= out_no[ea];
				}*/

			} // cell not empty
		} // tru cells in matrix

		//outputting
		output << time << ';' << rnarep::CellContent::no_replicators;
		double no;
		for(int ea = 0; ea <= par_noEA; ea++) {
			if(no = (double) out_no[ea]){
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
		- text file, tab separated, each line represents a grid point from 0th to last
			values:
			seq str mfe Pfold Pdeg no_sites R M type [activities] prev_type
		- rng binary state file
		*/
		std::string emptystring("N\tN\t0\t-1\t-1\t-1\t-1\t-1\t0");
		std::string filename;

		//prepare string for empty cells
		for(int ea = 0; ea < par_noEA; ea++) emptystring += "\t0";
		emptystring += "\t-1\n";

		//open outputs
		if(!savedir.length()) {
			std::cerr << "ERROR: No savedir inicialised! Please do run CompartPool::openOutputs() before saving!" << std::endl;
			return (1);
		}
		
		filename = savedir; 
		filename += '/'; 
		filename += std::to_string(time);
		filename += ".tsv";
		std::ofstream out(filename);
		
		if(!out){
			std::cerr << "ERROR: Could not open file for saving grid data!" << std::endl;
			return 2;
		}

		//going throught grid
		for(Compart *cell = matrix, *end = (Compart *) matrix + size ; cell != end; cell++){
			rnarep::CellContent *cellcont = cell->vals;
			//output values:
			// seq str mfe Pfold Pdeg no_sites R type [alphas]
			if(cellcont->empty) out << emptystring; 
			else {
				out 	<< *(cellcont->get_seq())
					<< '\t' << cellcont->get_str()
					<< '\t' << cellcont->get_mfe() 
					<< '\t' << cellcont->getPfold()
					<< '\t' << cellcont->Pdeg 
					<< '\t' << cellcont->get_no_sites()
					<< '\t' << cellcont->getR()
					<< '\t' << cell->M()
					<< '\t' << cellcont->get_type() ; 
				for(double *a = cellcont->geta(), *a_until = cellcont->geta() + par_noEA; a != a_until; a++){
					out << '\t' << *a;
				}

				out 	<< '\t' << cellcont->get_prev_type() << '\n';
			}

			//out.flush();
		}

		out.close();

		//saving rng state
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


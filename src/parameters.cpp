#include "parameters.h" 
#include "include/cm.h"

using namespace std;

//parameters harder to change
				int par_noEA=3;

//parameters to change with Args()
//ca params
				int par_quit = cmdv::QuitCond::qnone; // quitting condition. See: cdv::QuitCond
				int par_maxtime = 2000000; //length of simulation

				double par_claimNorep = 0.0; //claim of no replication - compared with M*R: it scales the possibility of a replicator daemon just hangin around not doing anyting
				double par_init_grid = 0.0; //
				int par_splitfrom = 50; //size when cells split
				int par_num_input_content = 25; //number of replicators read from file to cells
				int par_poolsize = 300; //number of replicators read from file to cells

//seed
				int par_seed = -1;
				int par_seed_plus = 0;
				char par_seed_file[255] = "\0";
//outputting
				int par_output_interval = 1000; // this often it writes averaged output
				int par_save_interval = 1000000; // thi often it saves all replicators to files

				char par_ID[255] = "test\0"; //name of simulation -> name of results folder
				char par_str_pool[255] = "IN/mapping.txt"; // single file containing decription of enzymatic activities
				char par_outdir[255] = "OUT"; // name of output folder
				char par_output_filename[255] = "output.csv"; // the name of file, where average output is written
				char par_savedir[255] = "SAVE"; //
				char par_load[255] = "\0"; // initial state
				char par_bubbles[255] = "\0";

//mutation rates
				double par_substitution = 0.005; //
				double par_insertion = 0.0005; //
				double par_deletion = 0.0005; //

//equations
/*
Pfold = exp(-cE) / (1 + exp(-cE))
e is mfe of replicator
*/
				double par_c = -0.3; //minus c!!!

/*
a = Pfold * alpha / m^sigma
alpha is activity of a motif
m is number of motifs in a replicator
*/
				double par_sigma = 1.1; //
				double par_gc_bonus = -0.3; //

/*
R_s = g * (l * 1 - P_{fold}) / (b1 + b2 * length)
length is length of the sequence
P_fold is calculated (see above)
*/
				double par_g = 10; //
				double par_b1 = 0.75; //
				double par_b2 = 0.005; //
				double par_ll = 2; // l + 1 in equation Rs 

/*
Pdeg = 0.9 - 0.8 * e / Emin
e is mfe of replicator
*/
				double par_Emin = -25.0; //

				double par_rangePdeg = 0.8; //
				double par_minPdeg = 0.1; //
				double par_flexPdeg = 0.2; //






//output parameters to file
int paramsToFile(const char* filename){
	// open a file in write mode.
	std::fstream paramfile(filename, std::fstream::out);
	
//	std::cout << "Printing parameters to file: " << filename << std::endl;	
	
	//write vienna rna version to file
	std::string comm("RNAfold --version >> ");
	comm += filename;
	system(comm.c_str());


	//paramfile << "RNAversion " << VRNA_VERSION << std::endl;
	paramfile << "MAXLEN " << MAXLEN << std::endl;
	paramfile << "par_noEA " << par_noEA << std::endl;
	paramfile << "par_quit " << par_quit << std::endl;
	paramfile << "par_maxtime " << par_maxtime << std::endl;
	paramfile << "par_poolsize " << par_poolsize << std::endl;
	paramfile << "par_splitfrom " << par_splitfrom << std::endl;
	paramfile << "par_claimNorep " << par_claimNorep << std::endl;
	paramfile << "par_num_input_content " << par_num_input_content << std::endl;
	paramfile << "par_output_interval " << par_output_interval << std::endl;
	paramfile << "par_save_interval " << par_save_interval << std::endl;
	paramfile << "par_seed " << par_seed << std::endl;
	paramfile << "par_seed_plus " << par_seed_plus << std::endl;
	paramfile << "par_ID " << par_ID << std::endl;
	paramfile << "par_str_pool " << par_str_pool << std::endl;
	paramfile << "par_outdir " << par_outdir << std::endl;
	paramfile << "par_output_filename " << par_output_filename << std::endl;
	paramfile << "par_savedir " << par_savedir << std::endl;
	paramfile << "par_load " << par_load << std::endl;
	paramfile << "par_seed_file " << par_seed_file << std::endl;
	paramfile << "par_bubbles " << par_bubbles << std::endl;
	//paramfile << "par_death " << par_death << std::endl;
	//paramfile << "par_diffusion_rate " << par_diffusion_rate << std::endl;
	paramfile << "par_init_grid " << par_init_grid << std::endl;
	paramfile << "par_ll " << par_ll << std::endl;
	paramfile << "par_sigma " << par_sigma << std::endl;
	//paramfile << "par_claimEmpty " << par_claimEmpty << std::endl;
	paramfile << "par_substitution " << par_substitution << std::endl;
	paramfile << "par_insertion " << par_insertion << std::endl;
	paramfile << "par_deletion " << par_deletion << std::endl;
	paramfile << "par_g " << par_g << std::endl;
	paramfile << "par_b1 " << par_b1 << std::endl;
	paramfile << "par_b2 " << par_b2 << std::endl;
	paramfile << "par_c " << par_c << std::endl;
	paramfile << "par_Emin " << par_Emin << std::endl;
	paramfile << "par_rangePdeg " << par_rangePdeg << std::endl;
	paramfile << "par_minPdeg " << par_minPdeg << std::endl;
	paramfile << "par_flexPdeg " << par_flexPdeg << std::endl;
	//paramfile << "par_Nmet " << par_Nmet << std::endl;
	//paramfile << "par_Nrep " << par_Nrep << std::endl;
	paramfile << "par_gc_bonus " << par_gc_bonus << std::endl;
	paramfile << "versioninfo " << versioninfo << std::endl;
	

	paramfile.close();
	return(0);
}


//parameter modifying function
int Args(int argc, char **argv)
{
  int i;
  char option;
  //  char temp[BUFSIZE];
  /* parameters may be overridden here. */
  for (i = 1; i < argc; i++) {
	if(argv[i][0] == '-'){
		option = argv[i][1];
		if( option == '-'){ //long expression
			//if(!strcmp(argv[i], "--par_death")) option = 'k';
			if(!strcmp(argv[i], "--par_diffusion_rate")) option = 'D';
			else if(!strcmp(argv[i], "--par_maxtime")) option = 'T';
			else if(!strcmp(argv[i], "--par_num_input_content")) option = 'n';
			else if(!strcmp(argv[i], "--par_poolsize")) option = 'N';
			else if(!strcmp(argv[i], "--par_splitfrom")) option = 'S';
			//else if(!strcmp(argv[i], "--par_nrow")) option = 'R';
			else if(!strcmp(argv[i], "--par_output_interval")) option = 'o';
			else if(!strcmp(argv[i], "--par_save_interval")) option = 'w';
			else if(!strcmp(argv[i], "--par_seed")) option = 'y';
			else if(!strcmp(argv[i], "--par_seed_plus")) option = '+';
			else if(!strcmp(argv[i], "--par_ID")) option = 'I';
			else if(!strcmp(argv[i], "--par_str_pool")) option = 'P';
			else if(!strcmp(argv[i], "--par_outdir")) option = 'O';
			else if(!strcmp(argv[i], "--par_output_filename")) option = 'F';
			else if(!strcmp(argv[i], "--par_savedir")) option = 'A';
			else if(!strcmp(argv[i], "--par_load")) option = 'L';
			else if(!strcmp(argv[i], "--par_seed_file")) option = 'f';
			//else if(!strcmp(argv[i], "--par_init_grid")) option = 'S';
			else if(!strcmp(argv[i], "--par_ll")) option = 'l';
			else if(!strcmp(argv[i], "--par_sigma")) option = 'G';
			else if(!strcmp(argv[i], "--par_claimNorep")) option = 'M';
			else if(!strcmp(argv[i], "--par_substitution")) option = 's';
			else if(!strcmp(argv[i], "--par_insertion")) option = 'i';
			else if(!strcmp(argv[i], "--par_deletion")) option = 'd';
			else if(!strcmp(argv[i], "--par_g")) option = 'g';
			else if(!strcmp(argv[i], "--par_b1")) option = '1';
			else if(!strcmp(argv[i], "--par_b2")) option = '2';
			else if(!strcmp(argv[i], "--par_c")) option = 'c';
			else if(!strcmp(argv[i], "--par_Emin")) option = 'e';
			//else if(!strcmp(argv[i], "--par_Nmet")) option = 'm';
			//else if(!strcmp(argv[i], "--par_Nrep")) option = 'r';
			else if(!strcmp(argv[i], "--par_gc_bonus")) option = 'U';
			else if(!strcmp(argv[i], "--par_rangePdeg")) option = 'x';
			else if(!strcmp(argv[i], "--par_minPdeg")) option = 'X';
			else if(!strcmp(argv[i], "--par_flexPdeg")) option = 'k';
			else if(!strcmp(argv[i], "--par_bubbles")) option = 'b';
			else if(!strcmp(argv[i], "--par_quit")) option = 'Q';
			else {
				std::cerr << "ERROR at reading argoments: not valid argoment (" << argv[i] << ")!" << std::endl;
				return -1;
			}
		}
		switch(option){
			// double
//			case 'r':
//				if (++i == argc) return 1;
//				par_Nrep = atof(argv[i]);
//				//if(par_Nrep > 0 ) {
//				//	std::cerr << "ERROR at reading argoments: option " << option << ": par_Nrep have to be ve!" << std::endl;
//				//	return(-1);
//				//}
//				continue;

//			case 'm':
//				if (++i == argc) return 1;
//				par_Nmet = atof(argv[i]);
//				//if(par_Nmet > 0 ) {
//				//	std::cerr << "ERROR at reading argoments: option " << option << ": par_Nmet have to be negative!" << std::endl;
//				//	return(-1);
//				//}
//				continue;

			case 'e':
				if (++i == argc) return 1;
				par_Emin = atof(argv[i]);
				if(par_Emin > 0 ) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_Emin have to be negative!" << std::endl;
					return(-1);
				}
				continue;

			case 'c':
				if (++i == argc) return 1;
				par_c = atof(argv[i]);
				if(par_c > 0 ) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_c have to be negative (it is the negative of scaling factor c)!" << std::endl;
					return(-1);
				}
				continue;

			case 'x':
				if (++i == argc) return 1;
				par_rangePdeg = atof(argv[i]);
				if(par_rangePdeg < 0 || par_rangePdeg > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;
			
			
			case 'X':
				if (++i == argc) return 1;
				par_minPdeg = atof(argv[i]);
				if(par_minPdeg < 0 || par_minPdeg > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'k':
				if (++i == argc) return 1;
				par_flexPdeg = atof(argv[i]);
				/*if(par_flexPdeg < 0 || par_flexPdeg > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": have to be between 0 and 1!" << std::endl;
					return(-1);
				}*/
				continue;

			case '2':
				if (++i == argc) return 1;
				par_b2 = atof(argv[i]);
				if(par_b2 < 0 ) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_b2 have to be positive!" << std::endl;
					return(-1);
				}
				continue;

			case '1':
				if (++i == argc) return 1;
				par_b1 = atof(argv[i]);
				if(par_b1 < 0 ) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_b1 have to be positive!" << std::endl;
					return(-1);
				}
				continue;

			case 'g':
				if (++i == argc) return 1;
				par_g = atof(argv[i]);
				if(par_g < 0 ) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_g have to be larger than 0!" << std::endl;
					return(-1);
				}
				continue;

			case 'd':
				if (++i == argc) return 1;
				par_deletion = atof(argv[i]);
				if(par_deletion < 0 || par_deletion > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_deletion have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;

			case 'i':
				if (++i == argc) return 1;
				par_insertion = atof(argv[i]);
				if(par_insertion < 0 || par_insertion > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_insertion have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;

			case 's':
				if (++i == argc) return 1;
				par_substitution = atof(argv[i]);
				if(par_substitution < 0 || par_substitution > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_substitution have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;

			case 'M':
				if (++i == argc) return 1;
				par_claimNorep = atof(argv[i]);
				if(par_claimNorep < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": has to be positive!" << std::endl;
					return(-1);
				}
				continue;

			/*case 'k':
				if (++i == argc) return 1;
				par_death = atof(argv[i]);
				if(par_death < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_death cant be negative!" << std::endl;
					return(-1);
				}
				continue;*/

			case 'l':
				if (++i == argc) return 1;
				par_ll = atof(argv[i]);
				if(par_ll < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_ll cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'G':
				if (++i == argc) return 1;
				par_sigma = atof(argv[i]);
				if(par_sigma <= 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_sigma have to be higher than 1!" << std::endl;
					return(-1);
				}
				continue;

//			case 'D':
//				if (++i == argc) return 1;
//				par_diffusion_rate = atof(argv[i]);
//				if(par_diffusion_rate < 0) {
//					std::cerr << "ERROR at reading argoments: option " << option << ": par_diffusion_rate cant be negative!" << std::endl;
//					return(-1);
//				}
//				continue;
			
//			case 'S':
//				if (++i == argc) return 1;
//				par_init_grid = atof(argv[i]);
//				if(par_init_grid < 0 || par_init_grid > 1) {
//					std::cerr << "ERROR at reading argoments: option " << option << ": have to be between 0 and 1!" << std::endl;
//					return(-1);
//				}
//				continue;
			
			case 'U':
				if (++i == argc) return 1;
				par_gc_bonus = atof(argv[i]);
				if(par_gc_bonus < 0 || par_gc_bonus > 1) {
					std::cerr << "ERROR at reading argoments: option " << option << ": have to be between 0 and 1!" << std::endl;
					return(-1);
				}
				continue;
			
			// int
			case 'Q':
				if (++i == argc) return 1;
				par_quit = atoi(argv[i]);
				if(par_quit < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			case '+':
				if (++i == argc) return 1;
				par_seed_plus = atoi(argv[i]);
				continue;
			case 'S':
				if (++i == argc) return 1;
				par_splitfrom = atoi(argv[i]);
				if(par_splitfrom <= 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			case 'y':
				if (++i == argc) return 1;
				par_seed = atoi(argv[i]);
				if(par_seed < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			case 'T':
				if (++i == argc) return 1;
				par_maxtime = atoi(argv[i]);
				continue;
			
			case 'n':
				if (++i == argc) return 1;
				par_num_input_content = atoi(argv[i]);
				if(par_num_input_content <= 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_num_input_content cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'N':
				if (++i == argc) return 1;
				par_poolsize = atoi(argv[i]);
				if(par_poolsize<= 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_poolsize cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'o':
				if (++i == argc) return 1;
				par_output_interval = atoi(argv[i]);
				if(par_output_interval < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;

			case 'w':
				if (++i == argc) return 1;
				par_save_interval = atoi(argv[i]);
				if(par_save_interval < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;


			// char
			case 'f':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) strcpy(par_seed_file, argv[i]);
				continue;
				
			case 'b':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) 
					strcpy(par_bubbles, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": should be more than 0 char long!" << std::endl;
					return -1;
				}
				continue;

			case 'L':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) 
					strcpy(par_load, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": should be more than 0 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'A':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) 
					strcpy(par_savedir, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": should be more than 0 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'F':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) 
					strcpy(par_output_filename, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": should be more than 0 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'O':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 0 ) 
					strcpy(par_outdir, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_outdir should be more than 0 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'I':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 1 ) 
					strcpy(par_ID, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": ID should be more than 1 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'P':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 1 ) 
					strcpy(par_str_pool, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": ID should be more than 1 char long!" << std::endl;
					return -1;
				}
				continue;
			
			default:
				std::cerr << "ERROR at reading argoments: not valid argoment (" << option << ")!" << std::endl;
				return -1;
		}
	}
  }
  return 0;
}

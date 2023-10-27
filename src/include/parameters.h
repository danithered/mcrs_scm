#ifndef _MYPARAMS_
#define _MYPARAMS_
 
#include <iostream>
#include <cstring>
#include <fstream>

//extern "C" {
//#include <ViennaRNA/vrna_config.h>
//}
 
#define MAXLEN 300

const char versioninfo[255] = "ver 3.1.5 new output\0";

extern int par_noEA;

extern int par_quit;
extern int par_maxtime;
extern int par_poolsize;
extern int par_output_interval;
extern int par_save_interval;
extern int par_seed;
extern int par_seed_plus;
extern double par_claimNorep;
extern int par_splitfrom;
extern int par_num_input_content;

extern char par_ID[255];
extern char par_str_pool[255];

extern char par_outdir[255];
extern char par_output_filename[255];
extern char par_savedir[255];
extern char par_load[255];
extern char par_seed_file[255];
extern char par_bubbles[255];
extern char par_str_pool_copy[255];

extern double par_init_grid;
extern double par_ll;
extern double par_sigma;
extern double par_substitution;
extern double par_insertion;
extern double par_deletion;
extern double par_g;
extern double par_b1;
extern double par_b2;
extern double par_c;
extern double par_Emin;
extern double par_gc_bonus;
extern double par_rangePdeg;
extern double par_minPdeg;
extern double par_flexPdeg;


//Functions 
int paramsToFile(const char* filename);
int Args(int argc, char **argv);

#endif

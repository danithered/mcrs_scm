#include <iostream>
#include <vector>
#include <string.h>
#include "randomgen.h"

using namespace std;

int find_new_pos(std::vector<int> &vec, int max){
	int hanyadik = gsl_ran_binomial(r, 0.5, max) , p;

	if( vec.empty() ) return(hanyadik);
	
	for (std::vector<int>::iterator  d = vec.begin(); (!vec.empty())  && d != vec.end(); ++d) {
		if(hanyadik == *d) {
			hanyadik = find_new_pos(vec, max);
			break;
		}
	}
	return(hanyadik);
}

char bases[] = "AGCU"; 

int give_random_str(int str_length = 20, double prop_paired = 0.4, double prop_seq = 0.2){

	int N = gsl_ran_binomial(r, 0.5, str_length);

	int n_p = gsl_ran_binomial(r, prop_paired, N/2);
	int no_paired = n_p*2;
	std::string str;
	int n_points = N - no_paired;

	std::vector<int> pos;
	std::vector<char> base;

	if(no_paired){
		str.assign(no_paired, '0');
//		std::cout << "N=" << N << ", no_paired=" << no_paired << std::endl;

		str[0] = '(';

		if(0.0 < gsl_rng_uniform(r)){ //not interloop
			//arranging (s and )s
			int n_c=1, n_rc=0;
			std::vector<std::string::iterator> where;

			for(std::string::iterator ch=str.begin() + 1; ch!=str.end(); ++ch){
				//std::cout << "n_p: " << n_p << ", n_c: " << n_c << ", n_rc: " << n_rc << std::endl;
				//std::cout << ( (double) n_p-n_c)/( (double) n_p-n_rc) << std::endl;
				if(gsl_rng_uniform(r) < ( (double) n_p-n_c)/( (double) n_p-n_rc) ){
					*ch = '(';
					n_c++;
				}
				else{
					*ch = ')';
					n_rc++;
				}
			}
			//adding points to it
			//first adding loop bases
			for(std::string::iterator ch=str.begin(); (ch+1) != str.end(); ++ch){
				if(*ch == '(' && *(ch+1) == ')' ) {
					where.push_back(ch+1);
//					std::cout << "where: " << (int) ch - str.begin() << std::endl;
				}
			}
			for(; n_points > 0 && where.size() > 0 ; n_points-- ) {
				int which = gsl_rng_uniform_int(r, where.size());
				str.insert(where[which ], '.');
				for(auto w = where.begin(); w != where.end(); ++w) if( *w > where[which] ) (*w)++;
				where.erase(where.begin() + which);
			}
			//adding extra bases
			for(; n_points > 0; n_points--){
				str.insert(str.begin() + gsl_rng_uniform_int(r, str.length()-1) + 1, '.') ;
			}


		} 
		else { // interloop
		}

	}
	else {
		str.assign(N, '.');
//		std::cout << "N=" << N << ", no_paired=" << no_paired << std::endl;
	}

	std::cout << str << std::endl;

	if(n_points = N - no_paired) {

		for(int n_clues = gsl_ran_binomial(r, prop_seq, n_points); n_clues--;){
			//int hanyadik = gsl_ran_binomial(r, 0.5, n_points) , p;
			//for (auto d = pos.begin(); d != NULL && d != pos.end(); ++d) 
			//for( p=0; hanyadik = (str[p] == '.'):(--hanyadik)?hanyadik ; p++) {};
			//pos = ;
			pos.push_back( find_new_pos(pos, n_points-1) );
			base.push_back( bases[gsl_rng_uniform_int(r, 4)] );
//			std::cout << n_clues << ": " << pos.back() << base.back() << std::endl;
		}
		for(auto po = pos.begin(); po != pos.end(); ++po){
			int it, point_found;
			for(it = 0, point_found = -1 ; point_found < *po && it != str.length(); it++) if(str[it] == '.') point_found++; 
			*po = it;
//			std::cout << *po << std::endl;
		}
	}

	for(int i = 0; i < pos.size(); i++) std::cout << pos[i] << " " << base[i] << std::endl;
	//its output starts from 1! Not 0!

	return 0;
}

int main(int argc, char *argv[]){
	//initialise rng
	time_t timer;
	randomszam_inic( time(&timer) , r);

	give_random_str();








	//close rng
	gsl_rng_free(r);
	return 0;
}



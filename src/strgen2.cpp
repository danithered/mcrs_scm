#include <iostream>
#include <fstream>
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

int main(int argc, char *argv[]){
	//initialise rng
	time_t timer;
	randomszam_inic( time(&timer) , r);

	int num = 0;
	system("mkdir -p OUT/rstr");

	//simple loop
	for(int looplength = 5; looplength < 19; looplength += 2){
		num++;

		//create dir
		std::string command = "mkdir -p OUT/rstr/";
		command += std::to_string(num);
		system(command.c_str());

		//write out str
		command.clear();
		command = "OUT/rstr/";
		command += to_string(num);
		command += "/str.txt";
		std::ofstream file(command.c_str(), std::ofstream::out);

		file << '(';
		for(int l = looplength; l--;) file << '.';
		file << ')' << std::endl;

		for(int l = (looplength+1)/2; l <= (looplength+5)/2; l++){
			file << l << ' ' << bases[gsl_rng_uniform_int(r, 4)] << std::endl;
		}
		file.close();

	}
	//interloop
	for(int looplength = 5; looplength < 19; looplength += 2){
		num++;

		//create dir
		std::string command = "mkdir -p OUT/rstr/";
		command += std::to_string(num);
		system(command.c_str());

		//get random sequence to the middle of the structure
		std::vector<char> randombases = {'N','N','N'};
		for(auto b = randombases.begin(); b != randombases.end(); b++) {
			*b = bases[gsl_rng_uniform_int(r, 4)];
//			std::cout << randombases[0] << randombases[1] << randombases[2] << std::endl;
		}

		//write out str
		command.clear();
		command = "OUT/rstr/";
		command += to_string(num);
		command += "/str.txt";
		std::ofstream file(command.c_str(), std::ofstream::out);

		// left
		file << '(';
		for(int l = looplength; l--;) file << '.';
		file << '(' << std::endl;

		auto rb = randombases.begin();
		for(int l = (looplength+1)/2; l <= (looplength+5)/2; l++){
			file << l << ' ' << *(rb++) << std::endl;
		}
	//	file << std::endl;

		//right
		file << ')';
		for(int l = looplength; l--;) file << '.';
		file << ')' << std::endl;

		rb = randombases.begin();
		for(int l = (looplength+1)/2; l <= (looplength+5)/2; l++){
			file << l << ' ' << *(rb++) << std::endl;
		}
		file.close();

	}






	//close rng
	gsl_rng_free(r);
	return 0;
}



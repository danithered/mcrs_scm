#ifndef _DV_TOOLS_
#define _DV_TOOLS_
#include <vector>

namespace dvtools {

	inline int Rmod(int a, int b) {
		int r = a%b; 
		return r >= 0 ? r : r + b;
	}

	inline double fracpart(double x){
		return x - (int) x;
	}

	int brokenStickVals(double *values, int noChoices, double sum, double random); 
	int brokenStickVals(std::vector<double> &values, double random); 

	class quickPosVals {
		public:
			double last;

			//Constructor
			quickPosVals(int _max, double (*_f)(int));
			quickPosVals();

			//Copy Constructor
			quickPosVals(const quickPosVals &obj);

			//Destructor
			~quickPosVals();

			//indexing
			double &operator[](int i){
				if(i <= max) last = vals[i];
				else last = f(i);
				return last;
			}

			//give a new function
			void setFunc(double (*_f)(int) );

			//calculate for the given range
			void setMax(int newmax);

			//get max
			int getMax();

		private:
			double *vals;
			int max;
			double (*f)(int);
			
	};

}

#endif

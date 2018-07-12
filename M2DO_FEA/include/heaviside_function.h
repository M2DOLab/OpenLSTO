#ifndef M2DO_FEA_HEAVISIDE_FUNCTION_H
#define M2DO_FEA_HEAVISIDE_FUNCTION_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <cassert>
#include <set>

#include <../../vendor/eigen3/Eigen/Dense>
#include <../../vendor/eigen3/Eigen/Sparse>

using namespace std ;
using namespace Eigen ;

namespace M2DO_FEA {

	class HeavisideFunction {
		
		private:
			
			//

		public:

			// Properties:
			
			// Delta = step of the Heaviside function.
			// Beta = shift of the 0 level boundary.
			double delta, beta ;

			// Methods:

			HeavisideFunction () ;
			HeavisideFunction (double delta, double beta) ;
			void print () ;

			double value (double x) ;
			double grad  (double x) ;

	} ;

}

#endif
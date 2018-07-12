#ifndef M2DO_FEA_NODE_H
#define M2DO_FEA_NODE_H

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

	class Node {
		
		private:
			//

		public:
			// Properties:
			int spacedim ;
			int id ; // Global number.
			vector<double> coordinates ;
			vector<int> dof ;

			// Nodal material property.
			double property;
			
			// Methods:
			Node (int spacedim) ;
			Node (int spacedim, int id, vector<double> coordinates) ;
			void Print () ;
			vector<double> ReturnFirstNCoordinates (int n) ; // Return the first n coordinates.

	} ;

}

#endif
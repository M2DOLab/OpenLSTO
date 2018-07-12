#ifndef M2DO_FEA_QUADRATURE_H
#define M2DO_FEA_QUADRATURE_H

using namespace std ;
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

namespace M2DO_FEA {

	class GaussianQuadrature {
		
		private:
			//

		public:

			// Properties:
			int spacedim, order ;
			vector<double> eta ;  // Gauss point locations.
			vector<double> w ;	   // Gaussian weights.

			// Methods:
			GaussianQuadrature () ;
			GaussianQuadrature (int spacedim, int order) ;
			void Print () ;
			vector<double> UpdateEtaCounter(vector<double> & eta_count) ;

	} ;

}

#endif
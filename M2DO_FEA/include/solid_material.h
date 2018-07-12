#ifndef M2DO_FEA_SOLID_MATERIAL_H
#define M2DO_FEA_SOLID_MATERIAL_H

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

	class SolidMaterial {
		
		private:
			
			//

		public:

			// Properties:
			int spacedim ;
			double E   ; 	// Young's modulus.
			double nu  ; 	// Poisson's ratio.
			double rho ; 	// Density.
			double alphaT ;	// Coefficient of thermal expansion.
			double h   ; 	// Thickness

			MatrixXd C, V ;
			
			// Methods:
			SolidMaterial (int spacedim, double E, double nu, double rho = 0.0, double alphaT = 0.0, double h = 1.00) ;
			void Print () ;

	} ;

}

#endif
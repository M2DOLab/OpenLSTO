#ifndef M2DO_FEA_BOUNDARY_CONDITIONS_H
#define M2DO_FEA_BOUNDARY_CONDITIONS_H

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

	class DirichletBoundaryConditions {

		private:
			//

		public:

			// Properties:
			vector<int> dof ; // Total imposed DOF's.
			vector<int> dof_zeros; // Only DOF's == 0.
			int mesh_n_dof ;
			vector<int> reduced_dof_to_dof_map ;
			vector<int> dof_to_reduced_dof_map ;

			// Amplitudes of the boundary conditions.
			vector<double> amplitude; // To accommodate zero and non-zero values.

			// Methods:
			DirichletBoundaryConditions () ;
			DirichletBoundaryConditions (vector<int> & dof_in, vector<double> & amplitude_in, int mesh_n_dof_in) ;
			void Print () ;
			void MapReducedDofToDof () ;

	} ;

	class PointValues {

		private:
			//

		public:

			// Properties:
			vector<int>    dof ;
			vector<double> values ;

			// Methods:
			PointValues (vector<int> & dofs_in, vector<double> & values_in) ;
			void Print () ;

	} ;

}

#endif
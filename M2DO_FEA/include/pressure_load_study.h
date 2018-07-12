#ifndef M2DO_FEA_PRESSURE_LOAD_STUDY_H
#define M2DO_FEA_PRESSURE_LOAD_STUDY_H

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

#include "mesh.h"
#include "boundary_conditions.h"
#include "solid_element.h"
#include "acoustic_element.h"
#include "ma38_solver.h"
#include "ma57_solver.h"

using namespace std ;
using namespace Eigen ;

namespace M2DO_FEA {

	class PressureLoadStudy {
		
		/*
			[K]{u} = {f}
		*/

		private:
			// Count number of dofs at current acoustic-structure configuration.

		public:

			// Properties:

			Mesh & mesh ;

			
			SparseMatrix<double> K ;

			// VectorXd g ;

			VectorXd f ; // Force vector with Neumann boundary conditions.
			VectorXd f_reduced ;

			VectorXd u ;
			VectorXd u_reduced ;

			DirichletBoundaryConditions dirichlet_boundary_conditions ;

			// Methods:
			PressureLoadStudy (Mesh & mesh) ;
			void Print () ;

			void AddBoundaryConditions (DirichletBoundaryConditions) ;

			void AssembleF (PointValues &, bool time_it) ;
			MatrixXd CouplingElement2D (vector<int>) ;
			void AssembleKPressureLoad (bool time_it) ;
			void AssembleKFPressureLoadForOptimisation (bool time_it, vector<int>, vector<int>, vector<double>, vector<int>, vector<double>) ;
			void SolveWithHSLMA38 (bool print_output, bool time_it) ;
			void SolveWithHSLMA38_NodeBasedDOFs (bool print_output, bool time_it) ; // Solves already reduced DOF's.
			void SolveWithHSLMA57 (bool use_metis, bool print_output, bool time_it) ;
			// void SolveWithHSLMA57_TimesTransposed_NodeBasedDOFs (bool use_metis, bool print_output, bool time_it, vector<int>, vector<int>, vector<double>) ;

			// For cases with manually input pressure (no fluid being solved).
			void AssembleKFForOptimisation_NoFluid (bool time_it, vector<int>, vector<int>, vector<double>, vector<int>, vector<double>) ;
			void AssembleKFForOptimisation_NoFluid_NodalProperty (bool time_it, vector<int>, vector<int>, vector<double>, vector<int>, vector<double>, double);
			void AssemblePressureLoad_NoFluid (bool time_it, double, double, vector<double>, int, int) ;
			// void AssemblePressureLoad_EquivalentNodal_NoFluid (bool time_it, double, vector<vector<double> >) ; backup Wang et al 2015
			void AssemblePressureLoad_EquivalentNodal_NoFluid (bool time_it, double, vector<vector<double> >, double) ;

			// compute element compliances (strain energy).
			vector<double> ComputeElementsCompliances (bool time_it);
			vector<double> ComputeElementsCompliances_NodeBasedDOFs (bool time_it);

			// Save pressure load (when mechanical only).
			void SaveLoadTXT();
	} ;

}

#endif
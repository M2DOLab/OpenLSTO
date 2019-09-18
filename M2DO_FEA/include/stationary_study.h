#ifndef M2DO_FEA_STATIONARY_STUDY_H
#define M2DO_FEA_STATIONARY_STUDY_H

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

using namespace std ;
using namespace Eigen ;

namespace M2DO_FEA {

	class StationaryStudy {

		/*
			[K]{u} = {f}
		*/

		private:
			//

		public:

			// Properties:

			Mesh & mesh ;

			SparseMatrix<double> K ;

			VectorXd g ;

			VectorXd f ; // Force vector with Neumann boundary conditions.
			VectorXd f_reduced ;
			VectorXd f_i ; // Force vector to store adjoint or body loads.
			VectorXd f_i_reduced ;

			VectorXd u ;
			VectorXd u_reduced ;
			VectorXd u_i ;
			VectorXd u_i_reduced ;

			DirichletBoundaryConditions dirichlet_boundary_conditions ;

			// Methods:
			StationaryStudy (Mesh & mesh) ;
			void Print () ;

			void AddBoundaryConditions (DirichletBoundaryConditions) ;

			/*
				K = (grad v, c * grad u)_omega
			*/

			void AssembleKWithAreaFractions (bool time_it) ;
			void AssembleKWithNodalProperties (bool time_it, double) ; // (input = minimum property value).
			void AssembleF (PointValues &, bool time_it) ;
			void AddThermalExpansionLoad (double DeltaT, bool time_it) ;
			void SolveWithHSLMA57 (bool use_metis, bool print_output, bool time_it) ; // Solves [K] * {u_reduced} = {f_reduced}.
			void SolveWithCG () ; // Solves [K] * {u_reduced} = {f_reduced} using CG
			void SolveWithHSLMA57_f_i (bool use_metis, bool print_output, bool time_it); // solves the adjoint problem
			void AssembleF_i (MatrixXd lambda_i, vector<int> dof, bool time_it); // assebble the adjoint load
			void SolveWithCG_f_i ();
			// \brief
			//		Asssmbles the thermal expansion load for a non uniform temperature distribution
			// \param
			//	 A vector of the non uniform temperature at the nodes
			void AddThermalExpansionLoadNonUniform (VectorXd non_unif_temp, bool time_it) ;

			//Poisson Study

			SparseMatrix<double> Kc ;

			VectorXd fc ; // Force vector
			VectorXd fc_reduced ;
			VectorXd fconv_reduced ;

			VectorXd uc ;
			VectorXd uc_reduced ;

			void AssembleFc (PointValues &, bool time_it) ;
			void AssembleKcWithAreaFractions (bool time_it) ;
			void AssembleKcKnaturalconvWithAreaFractions (bool time_it, double T_ref, vector<vector<double> > boundary_segments, double scale) ; // input: outside reference temp, boundary_Segments, scaleLSMtoFEA
			void SolvePoissonWithCG () ; // Solves [K] * {u_reduced} = {f_reduced} using CG
			void SolvePoissonConvectionWithCG () ;
			void SolvePoissonWithHSLMA57(bool use_metis, bool print_output, bool time_it) ; // Solves [K] * {u_reduced} = {f_reduced}

			void AssembleKcKconvWithAreaFractions (const std::vector<double> &, const MatrixXd & , const VectorXd &) ; // input: convection coeff, average velocity, pressure
			void SolveAsymmetricPoissonWithHSLMA57(bool use_metis, bool print_output, bool time_it) ; // Solves [K] * {u_reduced} = {f_reduced}

			void AddPressureFieldLoad (bool time_it, vector<vector<double> > pressure_segments, double scale) ; // Add pressure load


	} ;

}

#endif

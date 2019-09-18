#ifndef M2DO_FEA_SOLID_ELEMENT_H
#define M2DO_FEA_SOLID_ELEMENT_H

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

	/*
		Forward-declare the Mesh class. In order to compute the
		Jacobian matrix (and perhaps other things), each element
		needs access to it's node coordinates (not just the node_ids).
		I think the best way is to give every element a reference to
		the mesh in which it belongs; thus, Mesh needs to be declared
		prior to Element.
	*/

	class Mesh ;

	class SolidElement {

		private:

			//

		public:

			// Properties:
			int spacedim, order ;
			vector<int> node_ids, dof ;
			Mesh & mesh ; // The mesh to which this element belongs.
			int material_id ;
			double area_fraction ;
			vector<double> centroid ;
			double radius ;

			LinearShapeFunction linear_shape_function ;
			GaussianQuadrature  quadrature ;

			// Methods:
			SolidElement (int spacedim, int order, Mesh & mesh) ;
			void Print () ;

			MatrixXd J (vector<double> & eta) ;
			MatrixXd B (vector<double> & eta) ;
			MatrixXd B_int () ;
			MatrixXd B_axisymmetric (vector<double> & eta, double radius) ;	// B matrix for axisymmetric case
			MatrixXd K () ;
			MatrixXd K_axisymmetric (double radius, int material_id) ;	// K matrix for axisymmetric case
			MatrixXd K_NodalProperties (VectorXd &, double) ; // K computed with nodal material interpolated properties (input = nodal properties, minimum density or minimum area fraction).
			vector<MatrixXd> dKdz (double) ; // derivative of K computed with nodal material interpolated properties (input = minimum density or minimum area fraction)..
			MatrixXd M () ;
			MatrixXd FThermalExpansion (double) ;
			VectorXd FSelfWeight (double, VectorXd) ;

			//! \brief
			//		Computes thermal load on a solid element based on an input
			//		temperature vector
			//! \param deltaT_vector
			//		vector of temperatures at the nodes
			MatrixXd FThermalExpansion (VectorXd deltaT_vector) ;

			VectorXd NaturalToPhysicalCoordinates (vector<double> & eta) ;
			MatrixXd PhysicalGaussPoissCoordinates () ;

	} ;

}

#endif

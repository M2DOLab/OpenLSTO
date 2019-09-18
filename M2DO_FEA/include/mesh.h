#ifndef M2DO_FEA_MESH_H
#define M2DO_FEA_MESH_H

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

#include "quadrature.h"
#include "linear_shape_function.h"
#include "node.h"
#include "solid_element.h"
#include "solid_material.h"

using namespace std ;
using namespace Eigen ;

namespace M2DO_FEA {

	class Mesh {

		private:
			// Compute centroids.
			void ComputeCentroids () ; // When only solid elements are used.
			void ComputeCentroidsAcousticStructure () ; // When acoustic and solid elements are used.

			// compute node connectivity and element neighbouring.
			void nodeConnectivity();
			void elementNeighbours();

		public:
			// Properties:
			int spacedim ;
			bool is_structured ;

			// Nodes.
			vector<Node> nodes ;
			vector<vector<int> > node_connectivity;

			// Elements and materials.
			vector<SolidElement> solid_elements ;
			vector<SolidMaterial> solid_materials ;
			vector<vector<int> > element_neighbours;
			// Initial fluids (remain fluids).
			vector<int> initial_fluids ;
			// Design domain (defined in the main file).
			vector<int> design_domain;

			// Vector to identify element types, when different elements are being used.
			vector<int> element_type; // 0 = void; 1 = solid; 2 = fluid;

			// Methods:
			Mesh () ;
			Mesh (int spacedim) ;
			void Print () ;

			// Purely structural problems.
			// Mesh.
			void MeshSolidHyperRectangle (vector<int> nel, MatrixXd mesh_box, int element_order, bool time_it) ;
			// Assign DOF's.
			void AssignDof () ;

			int n_dof ;
			int n_entries () ;
			// int nc_entries();

			vector<int> GetNodesByCoordinates (vector<double> coord, vector<double> tol) ;
			vector<int> dof (int node_id) ;
			vector<int> dof (vector<int> node_ids) ;
			vector<int> dof (vector<int> node_ids, vector<int> components) ;

			// Write node coordinates and element connectivity.
			void WriteMeshTXT () ;
			// Write area fractions.
			void WriteAreaFractionsTXT (int, int, int) ;
			// Write element types.
			void WriteElementTypesTXT (int, int, int) ;
			// Read element neighbours.
			void ReadNeighboursTXT(int, int);

			// Write nodal material properties.
			void saveNodalPropertiesVTK(const unsigned int&, int, int) const;
			// Write general nodal values (input nodal values, vector<double>).
			void saveNodalValuesVTK(const unsigned int&, int, int, vector<double>) const;


	} ;


}

#endif

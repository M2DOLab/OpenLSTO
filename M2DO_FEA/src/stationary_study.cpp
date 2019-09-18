#include "mesh.h"
#include "boundary_conditions.h"
#include "stationary_study.h"

using namespace M2DO_FEA ;

StationaryStudy :: StationaryStudy (Mesh & mesh) : mesh (mesh) {

	//

}

void StationaryStudy :: Print () {

	cout << "Stationary Study" ;

}

void StationaryStudy :: AddBoundaryConditions (DirichletBoundaryConditions bc_in) {

	dirichlet_boundary_conditions = bc_in ;

}

void StationaryStudy :: SolveWithCG_f_i () {

    ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;

    cg.compute(K) ;
    cg.setTolerance(1e-6);
    cg.setMaxIterations(10000);

    u_i_reduced = cg.solve (f_i_reduced) ;

    u_i = VectorXd::Zero (mesh.n_dof) ;
    int n_dof_reduced = mesh.n_dof - dirichlet_boundary_conditions.dof_zeros.size() ;

    for (int i = 0 ; i < n_dof_reduced ; ++i) {
        u_i (dirichlet_boundary_conditions.reduced_dof_to_dof_map[i]) += u_i_reduced(i) ;
    }

}

void StationaryStudy :: AssembleF (PointValues & point_values, bool time_it) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling {f} from point values ... " << flush ;

	}

	int n_dof = mesh.n_dof ;
	f = VectorXd::Zero (n_dof) ;

	int n_dof_reduced = n_dof - dirichlet_boundary_conditions.dof_zeros.size() ;
	f_reduced = VectorXd::Zero(n_dof_reduced) ;

	int reduced_dof_i ;

	for (int i = 0 ; i < point_values.dof.size() ; ++i) {

		f (point_values.dof[i]) += point_values.values[i] ;

		reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ point_values.dof[i] ] ;

		if(reduced_dof_i >= 0) f_reduced (reduced_dof_i) += point_values.values[i] ;

	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

// void StationaryStudy :: AssembleFc (PointValues & point_values, bool time_it) {

// 	auto t_start = chrono::high_resolution_clock::now() ;

// 	if (time_it) {

// 		cout << "\nAssembling {f} from point values ... " << flush ;

// 	}

// 	int nc_dof = mesh.nc_dof ;
// 	fc = VectorXd::Zero (nc_dof) ;

// 	int nc_dof_reduced = nc_dof - dirichlet_boundary_conditions.dof_zeros.size() ;
// 	fc_reduced = VectorXd::Zero(nc_dof_reduced) ;

// 	int reduced_dof_i ;

// 	for (int i = 0 ; i < point_values.dof.size() ; ++i) {

// 		fc (point_values.dof[i]) += point_values.values[i] ;

// 		reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ point_values.dof[i] ] ;

// 		if(reduced_dof_i >= 0) fc_reduced (reduced_dof_i) += point_values.values[i] ;

// 	}

// 	auto t_end = chrono::high_resolution_clock::now() ;

// 	if (time_it) {

// 		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

// 	}

// }

void StationaryStudy :: AssembleKWithAreaFractions (bool time_it) {

	/*
		Here we build a reduced K matrix; the Dirichlet dof's are not included.
	*/

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling [K] using area fraction method ... " << flush ;

	}

	int n_dof = mesh.n_dof - dirichlet_boundary_conditions.dof_zeros.size() ;
	int reduced_dof_i, reduced_dof_j ;

	typedef Triplet<double> T ;
	vector<T> triplet_list ;

	/*
		Over-size the triplet list to avoid having to resize
		during the loop. The mesh.n_entries() gives an upper
		bound.
	*/

	triplet_list.reserve (mesh.n_entries()) ;

	MatrixXd K_e ;

	/*
		Solid elements:
	*/

	for (int k = 0 ; k < mesh.solid_elements.size() ; ++k) {

		auto && element = mesh.solid_elements[k] ;

		/*
			This gives the global dof numbers for the element:
		*/

		vector<int> dof = element.dof ;

		if ( k == 0 or not mesh.is_structured ) {

			K_e = element.K() ;

		}

		for (int i = 0 ; i < dof.size() ; ++i) {

			// Convert global to reduced dof:
			reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[i] ] ;

			if ( reduced_dof_i >= 0 ) {

				for (int j = 0 ; j < dof.size() ; ++j) {

					// Convert global to reduced dof:
					reduced_dof_j = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[j] ] ;

					if ( reduced_dof_j >= 0 ) {

						// Now push to the triplet:
						triplet_list.push_back( T(reduced_dof_i, reduced_dof_j, element.area_fraction * K_e (i, j)) ) ;

					}

				}

			}

		}

	} // for solid elements.

	K.resize (n_dof, n_dof) ;
	K.setFromTriplets (triplet_list.begin(), triplet_list.end()) ;

	/*
 		Impose non homogeneous boundary conditions.
 	*/

	// Penalty for imposing non homogeneous dirichlet condition.
	double penalty = 1e20;

	// Impose non homogeneous dirichlet conditions.
	for (int i = 0 ; i < dirichlet_boundary_conditions.dof.size() ; ++i)
	{
		if (dirichlet_boundary_conditions.amplitude[i] != 0)
		{
			// Current dof.
			int dofpen = dirichlet_boundary_conditions.dof[i];

			// Dof in the reduced_dof matrices.
			int reduced_dofpen = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dofpen ];

			// Changing to equivalent stiffness.
			K.coeffRef(reduced_dofpen,reduced_dofpen) = penalty;

			// Computing equivalent force.
			f_reduced(reduced_dofpen) += penalty*dirichlet_boundary_conditions.amplitude[i];
		}
	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void StationaryStudy :: AssembleKWithNodalProperties (bool time_it, double min_property) {

	/*
		Here we build a reduced K matrix; the Dirichlet dof's are not included.
	*/

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling [K] using area fraction method ... " << flush ;

	}

	int n_dof = mesh.n_dof - dirichlet_boundary_conditions.dof_zeros.size() ;
	int reduced_dof_i, reduced_dof_j ;

	typedef Triplet<double> T ;
	vector<T> triplet_list ;

	/*
		Over-size the triplet list to avoid having to resize
		during the loop. The mesh.n_entries() gives an upper
		bound.
	*/

	triplet_list.reserve (mesh.n_entries()) ;

	MatrixXd K_e ;

	/*
		Solid elements:
	*/

	for (int k = 0 ; k < mesh.solid_elements.size() ; ++k) {

		auto && element = mesh.solid_elements[k] ;

		/*
			Compute element stiffness matrix with nodal material interpolated properties.
		*/

		VectorXd nodal_properties = VectorXd::Zero (4,1);
		nodal_properties(0) = mesh.nodes[element.node_ids[0]].property;
		nodal_properties(1) = mesh.nodes[element.node_ids[1]].property;
		nodal_properties(2) = mesh.nodes[element.node_ids[2]].property;
		nodal_properties(3) = mesh.nodes[element.node_ids[3]].property;

		MatrixXd K_e = element.K_NodalProperties(nodal_properties, min_property) ;

		/*
			This gives the global dof numbers for the element:
		*/

		vector<int> dof = element.dof ;

		/*
			Build triplet for matrix:
		*/

		for (int i = 0 ; i < dof.size() ; ++i) {

			// Convert global to reduced dof:
			reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[i] ] ;

			if ( reduced_dof_i >= 0 ) {

				for (int j = 0 ; j < dof.size() ; ++j) {

					// Convert global to reduced dof:
					reduced_dof_j = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[j] ] ;

					if ( reduced_dof_j >= 0 ) {

						// Now push to the triplet:
						triplet_list.push_back( T(reduced_dof_i, reduced_dof_j, K_e (i, j)) ) ;

					}

				}

			}

		}

	} // for solid elements.

	K.resize (n_dof, n_dof) ;
	K.setFromTriplets (triplet_list.begin(), triplet_list.end()) ;

	/*
 		Impose non homogeneous boundary conditions.
 	*/

	// Penalty for imposing non homogeneous dirichlet condition.
	double penalty = 1e10;

	// Impose non homogeneous dirichlet conditions.
	for (int i = 0 ; i < dirichlet_boundary_conditions.dof.size() ; ++i)
	{
		if (dirichlet_boundary_conditions.amplitude[i] != 0)
		{
			// Current dof.
			int dofpen = dirichlet_boundary_conditions.dof[i];

			// Dof in the reduced_dof matrices.
			int reduced_dofpen = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dofpen ];

			// Computing equivalent force.
			f_reduced(reduced_dofpen) += penalty*K.coeffRef(reduced_dofpen,reduced_dofpen)*dirichlet_boundary_conditions.amplitude[i];

			// Changing to equivalent stiffness.
			K.coeffRef(reduced_dofpen,reduced_dofpen) = penalty*K.coeffRef(reduced_dofpen,reduced_dofpen);
		}
	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void StationaryStudy :: SolveWithCG () {

	ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;

	cg.compute(K) ;
	cg.setTolerance(1e-6);
	cg.setMaxIterations(10000);
	u_reduced = cg.solve(f_reduced) ;
	//std::cout << "# FEA Iterations: " << cg.iterations() << std::endl;

	u = VectorXd::Zero (mesh.n_dof) ;
	int n_dof_reduced = mesh.n_dof - dirichlet_boundary_conditions.dof.size() ;

	for (int i = 0 ; i < n_dof_reduced ; ++i) {
		u (dirichlet_boundary_conditions.reduced_dof_to_dof_map[i]) += u_reduced(i) ;
	}

}

void StationaryStudy :: AssembleF_i (MatrixXd lambda_i, vector<int> dof, bool time_it) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling {f} from point values ... " << flush ;

	}

	// Initializing dof_i
	int reduced_dof_i ;

	// Loop through elemento dofs size.
	for (int i = 0 ; i < dof.size() ; ++i) {

		// Assemble lambda_i in f_i.
		f_i(dof[i]) += lambda_i(i);

		// cout << dof[i] << "\t";
		// Check if dof is not constrained, then assemble f_i_reduced.
		if (find(dirichlet_boundary_conditions.dof_zeros.begin(), dirichlet_boundary_conditions.dof_zeros.end(), dof[i]) != dirichlet_boundary_conditions.dof_zeros.end())
		{
			// dof is contrained.
		}
		else
		{
			// cout << dof[i] << "\t";
			// Compute reduced dofs.
			reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[i] ];

			// Assemble lambda_i_reduced in f_i_reduced.
			f_i_reduced(reduced_dof_i) += lambda_i(i);
		}
		// cout << "\n";
	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void StationaryStudy :: AddPressureFieldLoad (bool time_it, vector<vector<double> > pressure_segments, double scale) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nAssembling {f} from point values ... " << flush ;

	}

	int reduced_dof_i ;

	// Mesh parameters.
	double elsize = 1.0; // Always one in this LSM.

	// Loop through all pressure segments.
	for (int i = 0; i < pressure_segments.size(); i++)
	{
		// Computing geometric parameters.
		double dx = (pressure_segments[i][2]-pressure_segments[i][0]);
		double dy = (pressure_segments[i][3]-pressure_segments[i][1]);
		double xmean = dx/2;
		double ymean = dy/2;
		MatrixXd Phi = MatrixXd::Zero (2, 1) ;

	    	// Use normal vector from level set.
	    	Phi << pressure_segments[i][6],
	    		   pressure_segments[i][7];

		// Jacobian
		double jacobian = sqrt(pow(xmean,2)+pow(ymean,2));
		
		// Computing modified Gauss points (using one point integration, gi = 0, wi = 2);
		vector<double> g_mod; g_mod.resize(2);

		// Qsi = (x-xc)/a. Hard-coded for gi = 0;
		g_mod[0] = (pressure_segments[i][2]+pressure_segments[i][0])/2;
		g_mod[0] -= mesh.solid_elements[int(pressure_segments[i][5])].centroid[0]*scale;
		g_mod[0] /= 0.5*elsize;
		// Eta = (y-yc)/b.
		g_mod[1] = (pressure_segments[i][3]+pressure_segments[i][1])/2;
		g_mod[1] -= mesh.solid_elements[int(pressure_segments[i][5])].centroid[1]*scale;
		g_mod[1] /= 0.5*elsize;

		// Shape functions.
		MatrixXd N_shape_A = MatrixXd::Zero (2, 8) ;
	    N_shape_A << 0.25*(1-g_mod[0])*(1-g_mod[1]), 0, 0.25*(1+g_mod[0])*(1-g_mod[1]), 0, 0.25*(1+g_mod[0])*(1+g_mod[1]), 0, 0.25*(1-g_mod[0])*(1+g_mod[1]), 0,
		             0, 0.25*(1-g_mod[0])*(1-g_mod[1]), 0, 0.25*(1+g_mod[0])*(1-g_mod[1]), 0, 0.25*(1+g_mod[0])*(1+g_mod[1]), 0, 0.25*(1-g_mod[0])*(1+g_mod[1]);

		// Directional (coupling) matrix. = wi*NA^T*Phi*J. Hard-coded for wi = 2;
		VectorXd Le = 2*N_shape_A.transpose()*Phi*jacobian;

		// Assemble pressure vector!!
		int reduced_dof;

		// Element DOF's.
		vector<int> dof = mesh.solid_elements[int(pressure_segments[i][5])].dof ;

		for (int j = 0 ; j < dof.size() ; ++j)
		{
			f (dof[j]) += Le(j)*pressure_segments[i][8]/scale ;

			if (dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[j] ] > 0)
			{
				reduced_dof_i = dirichlet_boundary_conditions.dof_to_reduced_dof_map[ dof[j] ] ;

				f_reduced (reduced_dof_i) += Le(j)*pressure_segments[i][8]/scale;
			}
		}
	}

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}
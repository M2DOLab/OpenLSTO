#include "M2DO_FEA.h"
#include "lsm_3d.h"
#include "lsm_opti_3d.h"

using namespace std ;
using namespace Eigen ;

namespace FEA = M2DO_FEA ;

int main () {

	/*
		==========================================================================
		Set up physics objects
		==========================================================================
	*/

	/*
		Dimensionality of problem:
	*/

	const int spacedim = 3 ;

	/*
		FEA & level set mesh parameters:
	*/

    const unsigned int nelx = 40, nely = 8, nelz = 40;

    /*
        Create an FEA mesh object.
    */

    FEA::Mesh fea_mesh (spacedim) ;

    /*
        Mesh a hyper rectangle.
    */

    MatrixXd fea_box (8, 3) ;

		fea_box <<   0.0, 0.0, 0.0,
	          nelx, 0.0 , 0.0,
	          nelx, nely , 0.0,
	          0.0, nely,0.0 ,
	          0.0, 0.0 ,nelz,
	          nelx, 0.0 ,nelz,
	          nelx, nely ,nelz,
	          0.0, nely,nelz;

    vector<int> nel = {nelx, nely , nelz} ;

    int element_order = 3 ;
    fea_mesh.MeshSolidHyperRectangle (nel, fea_box, element_order, false) ;
    fea_mesh.is_structured = true ;
    fea_mesh.AssignDof () ;

	/*
		Add material properties:
	*/

    fea_mesh.solid_materials.push_back (FEA::SolidMaterial (spacedim, 1.0, 0.3, 1.0)) ;

    /*
        Next we specify that we will undertake a stationary study, which takes the form [K]{u} = {f}.
    */

    FEA::StationaryStudy fea_study (fea_mesh) ;

	/*
		Add a homogeneous Dirichlet boundary condition (fix some nodes).
	*/

    vector<int> fixed_nodes = fea_mesh.GetNodesByCoordinates ({0.0, 0.0 , nelz}, {1e9, 1e9, 1e-12}) ;
    vector<int> fixed_dof = fea_mesh.dof (fixed_nodes) ;

    vector<double> amplitude (fixed_dof.size(),0.0) ; // Values equal to zero.

    fea_study.AddBoundaryConditions (FEA::DirichletBoundaryConditions (fixed_dof, amplitude, fea_mesh.n_dof)) ;

    /*
        Apply load.
    */

    vector<int> load_node = fea_mesh.GetNodesByCoordinates ({nelx, 0.5*nely,0.4*nelz}, {2.1, 2.1, 1.e-12}) ;
    vector<int> load_dof  = fea_mesh.dof (load_node) ;


    vector<double> load_val (load_node.size() * spacedim) ;

    for (int i = 0 ; i < load_node.size() ; ++i) {

        load_val[spacedim*i]   = 0.0 ;
        load_val[spacedim*i+1] = 0.0 ;
				load_val[spacedim*i+2] = -1.0 ;

    }

    FEA::PointValues point_load (load_dof, load_val) ;
    fea_study.AssembleF (point_load, false) ;

		// Create sensitivity analysis instance.
    FEA::SensitivityAnalysis sens(fea_study) ;

		// sensitivity properties
		double p_norm = 16.0;
		int sens_type = 1;
		double least_squares_radius = 2.0;

		/*
			==========================================================================
			Set Level set properties
			==========================================================================
		*/

		//  Create an object of the levelset class
		LevelSet3D level_set_3d;

		//Declare box dimensions and initialize the box
		std::vector<double> LS2FEmap(3,1);

		uint box_x = nelx*LS2FEmap[0];
		uint box_y = nely*LS2FEmap[1];
		uint box_z = nelz*LS2FEmap[2];

		// dimensions of the branches of the L beam
		level_set_3d.nx1 = 0.4*nelx;
		level_set_3d.nz1 = 0.4*nelz;

		level_set_3d.SetBoxDimensions(box_x,box_y,box_z);// Set up dimensions

		//Gotta define the pointers to phi and grid_vel in the main, otherwise it gets deleted!
		level_set_3d.phi = new mp4Vector[level_set_3d.num_grid_pts];

  	level_set_3d.MakeLBeam(); // Initialize L beam

		level_set_3d.num_gauss_pts = 2;

		/*
			==========================================================================
			Set sensitivity data properties
			==========================================================================
		*/

		// Create a sensitivity object
		SensitivityData SensData;

		double MaxVol = 64.0*0.5; // in percentage
	  SensData.MaxVol = MaxVol;

	  std::vector<double> UB(2,0);
	  std::vector<double> LB(2,0);

		// pass mesh dimensions to sensitivity data
	  SensData.nx = box_x;
	  SensData.ny = box_y;
	  SensData.nz = box_z;

		SensData.nx1 = level_set_3d.nx1;
	  SensData.nz1 = level_set_3d.nz1;

	  SensData.LS2FEmap = LS2FEmap;

	  SensData.LB = LB;
	  SensData.UB = UB;

	  double move_limit = 0.2;
	  SensData.move_limit = move_limit;

		/*
			==========================================================================
			Begin iterations here
			==========================================================================
		*/

		int n_iterations = 0;
		int max_iterations = 100;

    while (n_iterations < max_iterations) {

    	++n_iterations ;

			// Discretize boundary using Marching Cubes
			if(n_iterations > 1) delete[] level_set_3d.triangle_array;
			level_set_3d.MarchingCubesWrapper();

			// Pass info to SensData
			SensData.bpointsize = level_set_3d.num_boundary_pts;
			SensData.bPoints = level_set_3d.boundary_pts_one_vector;
			SensData.pointAreas = level_set_3d.boundary_areas;
			SensData.bpointsize = level_set_3d.num_triangles;
			SensData.iter = n_iterations;

			// Set-up narrow band
			level_set_3d.SetupNarrowBand();


		 // Calculate volume fractions
			level_set_3d.CalculateVolumeFractions();

			SensData.volumeFractions = level_set_3d.volumefraction_vector;


			// -----------------------Assemble and Solve FEA

        // Assign area fractions.
        for (unsigned int i=0 ; i< fea_mesh.solid_elements.size() ; i++)
        {
            if (SensData.volumeFractions[i] < 1e-6) fea_mesh.solid_elements[i].area_fraction = 1e-6 ;
            else fea_mesh.solid_elements[i].area_fraction = SensData.volumeFractions[i] ;
        }

        /*
        Assemble stiffness matrix [K] using area fraction method:
        */

        fea_study.AssembleKWithAreaFractions (false) ;

        /*
            Solve equation:
        */

        fea_study.SolveWithCG () ;
				//fea_study.SolveWithHSLMA57 (false,false,false) ;

				SensData.compliance = fea_study.u.dot(fea_study.f);

				double compliance = SensData.compliance;

        // Compute stress sensitivities (stress*strain) at the Gauss points.
				sens.ComputeStressSensitivities3D (false,p_norm) ;

				SensData.bsens.resize(SensData.bpointsize);
				SensData.vsens.resize(SensData.bpointsize);

				for (int i=0 ; i < SensData.bpointsize ; i++)
        {
					// current boundary point
					std::vector<double> boundary_point(3);
	        boundary_point[0] = SensData.bPoints[3*i];
	        boundary_point[1] = SensData.bPoints[3*i+1];
	        boundary_point[2] = SensData.bPoints[3*i+2];

					// compute boundary sensitivity at this point

					sens.ComputeBoundarySensitivities(boundary_point, least_squares_radius , sens_type, p_norm) ;

					// Assign sensitivities.
					SensData.bsens[i] = sens.boundary_sensitivities[2*i] ;
					SensData.vsens[i] = -1 ;

					//symmetry
					 boundary_point[1] = nely - SensData.bPoints[3*i+1];
					 sens.ComputeBoundarySensitivities(boundary_point, least_squares_radius , sens_type, p_norm) ;

 					// Assign sensitivities.
 					SensData.bsens[i] += sens.boundary_sensitivities[2*i+1] ;
				}

				sens.boundary_sensitivities.clear();

				// Optimize boundary point movement
				PerformOptimization_Stress_LBeam(SensData);

				// Print stats
				cout<< "Iter: "<< n_iterations << " Obj: " << sens.objective << " MaxStress: " << sens.von_mises_max << " Constr: " << SensData.Vol/double(nelx *nely * nelz) << endl;


		    // Resize optimum velocities
		    level_set_3d.opt_vel.resize(level_set_3d.num_boundary_pts);


	  	  level_set_3d.opt_vel = SensData.opt_vel_nlopt;

				for (int i=0 ; i < SensData.bpointsize ; i++)
        {
					// current boundary point
					std::vector<double> boundary_point(3);
	        boundary_point[0] = SensData.bPoints[3*i];
	        boundary_point[1] = SensData.bPoints[3*i+1];
	        boundary_point[2] = SensData.bPoints[3*i+2];

					// assign positive values to points on the boundary
					if(boundary_point[0] >= nelx - 2 && boundary_point[2] >= 0.4*nelz - 2
					&& boundary_point[2] <= 0.4*nelz + 0.1) level_set_3d.opt_vel[i] = 0.2 ;

				}



			  // Extrapolate velocities
			  level_set_3d.ExtrapolateVelocities();

			 //Extend velocities using fast marching method
			 // FMM inside...
			  level_set_3d.indices_considered = level_set_3d.indices_considered_inside;
			  level_set_3d.FastMarchingMethod();

		    //FMM outside...
			  for(int i = 0; i < level_set_3d.num_grid_pts; i++)   level_set_3d.phi_temp[i] = -level_set_3d.phi_temp[i]; // flip sign
			  level_set_3d.indices_considered = level_set_3d.indices_considered_outside;
			  level_set_3d.FastMarchingMethod();
			  for(int i = 0; i < level_set_3d.num_grid_pts; i++) level_set_3d.phi_temp[i] = -level_set_3d.phi_temp[i]; // flip sign

				// print out stl file after every 2 iterations
				if(fmod( n_iterations , 2 ) == 1)
				{
					int box_smooth = 1;
					level_set_3d.WriteSTL(box_smooth);
			 	}


			 // Advect
			 level_set_3d.Advect_LBeam();




    }

		// write to stl
		int box_smooth = 1;
		level_set_3d.WriteSTL(box_smooth);

	/*
		Clean up
	*/
	delete[] level_set_3d.triangle_array;
	delete[] level_set_3d.phi;

	cout << "\nProgram complete.\n\n" ;

	return 0 ;
}

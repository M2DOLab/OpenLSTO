#include "M2DO_FEA.h"
#include "M2DO_LSM.h"

#include <sys/time.h>
#include <fstream>

using namespace std ;

namespace FEA = M2DO_FEA ;
namespace LSM = M2DO_LSM ;

#include "lsm_hole_insertion.hpp"
#include "boundary_hole.hpp"

double get_wall_time () {

	struct timeval time ;
	if (gettimeofday (&time,NULL)) {
		//  Handle error
		return 0 ;
	}

	return (double)time.tv_sec + (double)time.tv_usec * .000001 ;

}

double get_cpu_time(){

	return (double)clock() / CLOCKS_PER_SEC ;

}

int main () {

	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                SETTINGS FOR THE FINITE ELEMENT ANALYSIS                 //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		FEA Mesh:
	*/

	// FEA mesh object for 2D analysis:
	FEA::Mesh fea_mesh (2) ;

	// Number of elements in x and y directions:
	const unsigned int nelx = 160, nely = 80 ;

	// fea_box contains the (x,y) coordinates of 4 corner points of rectangle containing the mesh:
	MatrixXd fea_box(4,2);

	fea_box <<   0.0,   0.0,
					nelx,    0.0,
					nelx, nely,
					   0.0, nely;


	// Element Gauss integration order:
  	int element_order = 2 ;

  	// Create structured mesh and assign degrees of freedom:
	fea_mesh.MeshSolidHyperRectangle ({nelx, nely}, fea_box, element_order, false) ;
	fea_mesh.is_structured = true ;
	fea_mesh.AssignDof () ;

	/*
		Define material properties and add to mesh:
	*/

	double   E = 1.0 ; // Young's Modulus
	double  nu = 0.3 ; // Poisson's ratio
	double rho = 1.0 ; // Density

	fea_mesh.solid_materials.push_back (FEA::SolidMaterial (2, E, nu, rho)) ;

	/*
		Next we specify that we will undertake a stationary study, which takes the
		form [K]{u} = {f}:
	*/

 	FEA::StationaryStudy fea_study (fea_mesh) ;

	/*
		Define homogeneous Dirichlet boundary condition (fixed nodes) and add to study:
	*/

	// Example 1: cantilever beam

	// Select dof using a box centered at coord of size tol:
	vector<double>    coord = {0.0, 0.0}, tol = {1e-12, 1e10} ;
	vector<int> fixed_nodes = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	vector<int>   fixed_dof = fea_mesh.dof (fixed_nodes) ;

 	// Example 2: half of simply supported beam or MBB beam

 	// Left boundary condition
	// vector<double>        coord_left = {0.0, 0.0}, tol_left = {1e-12, 1e10} ;
	// vector<int>     fixed_nodes_left = fea_mesh.GetNodesByCoordinates (coord_left, tol_left) ;
	// vector<int> fixed_condition_left = {0} ; // set fixed in only the x direction.
	// vector<int>       fixed_dof_left = fea_mesh.dof (fixed_nodes_left, fixed_condition_left) ;

	// // Right boundary condition
	// vector<double>        coord_right = {nelx, 0.0}, tol_right = {1e-12, 1e-12} ;
	// vector<int>     fixed_nodes_right = fea_mesh.GetNodesByCoordinates(coord_right, tol_right) ;
	// vector<int> fixed_condition_right = {1} ; // set fixed in only the y direction.
	// vector<int>       fixed_dof_right = fea_mesh.dof(fixed_nodes_right, fixed_condition_right) ;

	// // Combine dofs into a single vector
	// vector<int> fixed_dof ;
	// fixed_dof.reserve(fixed_dof_left.size() + fixed_dof_right.size()) ;
	// fixed_dof.insert(fixed_dof.end(), fixed_dof_left.begin(), fixed_dof_left.end()) ;
	// fixed_dof.insert(fixed_dof.end(), fixed_dof_right.begin(), fixed_dof_right.end()) ;

	// Add boundary conditions to study:
	// fea_study.AddBoundaryConditions (FEA::DirichletBoundaryConditions (fixed_dof, fea_mesh.n_dof)) ;
	vector<double> amplitude (fixed_dof.size(),0.0) ; // Values equal to zero.

    fea_study.AddBoundaryConditions (FEA::DirichletBoundaryConditions (fixed_dof, amplitude, fea_mesh.n_dof)) ;

	/*
		Define a point load of (0, -0.5) at the point (nelx, 0.5*nely) and add to study:
	*/

	// Example 1: cantilever beam

	// Select dof using a box centered at coord of size tol:
	coord = {1.0*nelx, 0.5*nely}, tol = {1e-12, 1e-12} ;
	vector<int> load_node = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	vector<int>  load_dof = fea_mesh.dof (load_node) ;

	vector<double> load_val (load_node.size() * 2) ;
	for (int i = 0 ; i < load_node.size() ; ++i) {
		load_val[2*i]   = 0.00 ; // load component in x direction.
		load_val[2*i+1] = -0.5 ; // load component in y direction.
	}

	// Example 2: half of simply supported beam or MBB beam

	// coord = {0.0, nely}, tol = {1e-12, 1e-12} ;
	// vector<int>      load_node = fea_mesh.GetNodesByCoordinates (coord, tol) ;
	// vector<int> load_condition = {1} ; // apply load in only the y direction.
	// vector<int>       load_dof = fea_mesh.dof (load_node, load_condition) ;

	// vector<double> load_val (load_node.size()) ;
	// for (int i = 0 ; i < load_node.size() ; ++i) {
	// 	load_val[i] = -10.0; //load component in y direction
	// }

	// Add point load to study and assemble load vector {f}:
	FEA::PointValues point_load (load_dof, load_val) ;
	fea_study.AssembleF (point_load, false) ;

	/*
		FEA Solver:
	*/

	// Initialise guess solution for CG:
	vector<double> u_guess (fea_mesh.n_dof, 0.0) ;

	// Convergence tolerance:
	double cg_tolerence = 1.0e-6 ;

	// END OF SETTINGS FOR THE FINITE ELEMENT ANALYSIS.


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                  SETTINGS FOR THE SENSITIVITY ANALYSIS                  //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	FEA::SensitivityAnalysis sens (fea_study) ;

	// END OF SETTINGS FOR THE SENSITIVITY ANALYSIS


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                   SETTINGS FOR THE LEVEL SET METHOD                     //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Define LSM parameters:
	*/

	double    move_limit = 0.5 ;   // Maximum displacement per iteration in units of the mesh spacing.
	double    band_width = 6 ;     // Width of the narrow band.
	bool is_fixed_domain = false ; // Whether or not the domain boundary is fixed.

	/*
		Seed initial holes:
		In this example, we create five horizontal rows, each row alternating between
		four and five equally spaced holes, all of radius 5 units.
	*/

	vector<LSM::Hole> holes ;

	// First row with five holes:
	// holes.push_back (LSM::Hole (16, 14, 5)) ;
	// holes.push_back (LSM::Hole (48, 14, 5)) ;
	// holes.push_back (LSM::Hole (80, 14, 5)) ;
	// holes.push_back (LSM::Hole (112, 14, 5)) ;
	// holes.push_back (LSM::Hole (144, 14, 5)) ;

	// // Second row with four holes:
	// holes.push_back (LSM::Hole (32, 27, 5)) ;
	// holes.push_back (LSM::Hole (64, 27, 5)) ;
	// holes.push_back (LSM::Hole (96, 27, 5)) ;
	// holes.push_back (LSM::Hole (128, 27, 5)) ;

	// // Third row with five holes:
	// holes.push_back (LSM::Hole (16, 40, 5)) ;
	// holes.push_back (LSM::Hole (48, 40, 5)) ;
	// holes.push_back (LSM::Hole (80, 40, 5)) ;
	// holes.push_back (LSM::Hole (112, 40, 5)) ;
	// holes.push_back (LSM::Hole (144, 40, 5)) ;

	// // Fourth row with four holes:
	// holes.push_back (LSM::Hole (32, 53, 5)) ;
	// holes.push_back (LSM::Hole (64, 53, 5)) ;
	// holes.push_back (LSM::Hole (96, 53, 5)) ;
	// holes.push_back (LSM::Hole (128, 53, 5)) ;

	// // Fifth row with five holes:
	// holes.push_back (LSM::Hole (16, 66, 5)) ;
	// holes.push_back (LSM::Hole (48, 66, 5)) ;
	// holes.push_back (LSM::Hole (80, 66, 5)) ;
	// holes.push_back (LSM::Hole (112, 66, 5)) ;
	// holes.push_back (LSM::Hole (144, 66, 5)) ;

	// END OF SETTINGS FOR THE LEVEL SET METHOD


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                      SETTINGS FOR THE OPTIMIZATION                      //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Define parameters needed for optimization loop:
	*/

	double       max_time = 6000 ;   // maximum running time.
	int    max_iterations = 300 ;    // maximum number of iterations.
	double       max_area = 0.5 ;    // maximum material area.
	double       max_diff = 0.0001 ; // relative difference between iterations must be less than this value to reach convergence.

	/*
		Lambda values for the optimiser:
		These are reused, i.e. the solution from the current iteration is
		used as an estimate for the next, hence we declare the vector
		outside of the main loop.
	*/

	vector<double> lambdas (2) ;

	// END OF SETTINGS FOR THE OPTIMIZATION


	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                   LEVEL SET TOPOLOGY OPTIMIZATION LOOP                  //
	//               SHOULD NOT NEED TO EDIT FILE BEYOND THIS POINT            //
	//                                                                         //
	/////////////////////////////////////////////////////////////////////////////

	/*
		Create level set
	*/

	// Initialise the level set mesh (same resolution as the FE mesh):
	LSM::Mesh lsm_mesh (nelx, nely, false) ;

	double mesh_area = lsm_mesh.width * lsm_mesh.height ;

	// Initialise the level set object (from the hole vector):
	LSM::LevelSet level_set (lsm_mesh, holes, move_limit, band_width, is_fixed_domain) ;

	// Reinitialise the level set to a signed distance function:
	level_set.reinitialise () ;

	// Initialise the boundary object :
	LSM::Boundary boundary (level_set) ;

	// Initialise random number generator:
	LSM::MersenneTwister rng ;

	/*
		Optimization:
	*/

	// Declare parameters that will change within the optimization loop:
	unsigned int   n_reinit = 0 ;              // num cycles since signed dist reinitialisation.
	double         time     = 0 ;              // running time.
	vector<double> times, compliances, areas ; // time, compliance and area measurements.
	int            n_iterations = 0 ;          // iteration counter
	vector<double> objective_values ;          // vector to save objective history
	double         relative_difference = 1.0 ; // convergence criteria variable,

	// Initialise io object:
	LSM::InputOutput io ;


	// ----------------------------------------------------------
    // Define/Initialise parameters for hole insertion subroutine
    vector <double> h_index(lsm_mesh.nNodes);
    vector <double> h_lsf(lsm_mesh.nNodes);
    vector <bool> h_elem(lsm_mesh.nElements);
    vector <M2DO_LSM::h_Node> h_Nsens(lsm_mesh.nNodes);
    vector <M2DO_LSM::h_Node> h_Nsens_Temp(lsm_mesh.nNodes);
    vector <M2DO_LSM::h_Element> h_Esens(lsm_mesh.nElements);
    int h_count = 0;
    double holeCFL = 0.15;
    double lBand = 2.0;
    double h; // Size of level set element max (width, height)
    double h_bar;  // artifitial height for setting secondary level set function
    bool h_flag = false;  //
    bool isHole = true; // Is hole inseration function active?
    double newHoleAreaLimit = 0.02;
    //
    // Assign desired artificial height
    //
    h = ( lsm_mesh.width/nelx > lsm_mesh.height/nely ) ? (lsm_mesh.width/nelx) : (lsm_mesh.height/nely);
    h_bar = h;

    // LSM::Boundary_hole boundary_hole(levelSet) ;
    LSM::Boundary boundary_hole(level_set) ;
    // vector<double> gammas(2) ;
    //------------------------------------------------------------

	cout << "\nStarting compliance minimisation demo...\n\n" ;

	// Print output header:
	printf ("--------------------------------\n") ;
	printf ("%8s %12s %10s\n", "Iteration", "Compliance", "Area") ;
	printf ("--------------------------------\n") ;

	// Create a directory name for output
	system("mkdir -p results/history");
	system("mkdir -p results/level_set");
	system("mkdir -p results/area_fractions");
	system("mkdir -p results/boundary_segments");

	system("rm -f results/history/*.vtk results/history/*.txt");
	system("rm -f results/level_set/*.vtk results/level_set/*.txt");
	system("rm -f results/area_fractions/*.vtk results/area_fractions/*.txt");
	system("rm -f results/boundary_segments/*.vtk results/boundary_segments/*.txt");

	// Setup text file:
	ofstream history_file ;
	history_file.open ("results/history/history.txt", ios_base::app) ;
	history_file << "Iteration\tCompliance\tArea\n" ;
	history_file.close () ;

	while (n_iterations < max_iterations) {

		++n_iterations ;

		// Perform boundary discretisation:
		boundary.discretise (false, lambdas.size()) ;

		// Compute element area fractions:
		boundary.computeAreaFractions () ;

		// Assign area fractions:
		for (int i = 0 ; i < fea_mesh.solid_elements.size() ; i++) {

			if (lsm_mesh.elements[i].area < 1e-3) {
				fea_mesh.solid_elements[i].area_fraction = 1e-3 ;
			}

			else {
				fea_mesh.solid_elements[i].area_fraction = lsm_mesh.elements[i].area ;
			}

		}

		// Assemble stiffness matrix [K] using area fraction method:
		fea_study.AssembleKWithAreaFractions (false) ;


		// Solve equation using conjugant gradient (cg) method:
		fea_study.SolveWithCG();

		// Compute compliance sensitivities (stress*strain) at the Gauss points:
		sens.ComputeComplianceSensitivities (false) ;

		for (int i = 0 ; i < boundary.points.size() ; i++) {

			vector<double> boundary_point (2, 0.0) ;
			boundary_point[0] = boundary.points[i].coord.x ;
			boundary_point[1] = boundary.points[i].coord.y ;

			// Interpolate Gauss point sensitivities by least squares
			sens.ComputeBoundarySensitivities (boundary_point);

			// Assign sensitivities
			boundary.points[i].sensitivities[0] = -sens.boundary_sensitivities[i];
			boundary.points[i].sensitivities[1] = -1;

		}

		// clearing sens.boundarysens vector
		sens.boundary_sensitivities.clear () ;

		// Time step associated with the iteration
		double time_step ;

		// Constraint distance vector
		vector<double> constraint_distances ;

		// Push current distance from constraint violation into vector
		constraint_distances.push_back (mesh_area * max_area - boundary.area) ;

		/* Initialise the optimisation object

		The Optimise class is a lightweight object so there is no cost for
		reinitialising at every iteration. A smart compiler will optimise
		this anyway, i.e. the same memory space will be reused. It is better
		to place objects in the correct scope in order to aid readability
		and to avoid unintended name clashes, etc.
		*/

		LSM::Optimise optimise (boundary.points, time_step, move_limit) ;
		// LSM::Optimise optimise(boundary.points, constraint_distances,
		// 						lambdas, time_step, move_limit, false, true);



		// set up required parameters
		optimise.length_x = lsm_mesh.width ;
		optimise.length_y = lsm_mesh.height ;
		optimise.boundary_area = boundary.area ; // area of structure
		optimise.mesh_area = mesh_area ; // area of the entire mesh
		optimise.max_area = max_area ; // maximum area, i.e. area constraint

		// Perform the optimisation
		optimise.Solve_With_NewtonRaphson () ;

		optimise.get_lambdas(lambdas);
		// lambdas[0] = 0.1 * lambdas[0];
		// lambdas[1] = 0.1 * lambdas[1];

	    ///////////////////////////////////////////////////////////////////////////


	    /////////////////////////////////////////////////////////////////////////////
        //                                                                         //
        //       Conduct hole insertion scheme                                     //
        //         main procedures:                                                //
        //           1. Initialise nodes with a value to represent their           //
        //                artificial height, h_bar, in terms of element length, h  //
        //           2. Identify area/nodes that new holes can be inserted         //
        //           3. Calculate nodal sensitivities using least square method    //
        //           4. Update the value of the secondary level set function for   //
        //                hole insertable nodes                                    //
        //           5. Check whether new hole should be inserted if corresponding //
        //                node's value of the secondary level set function is less //
        //                than 0. If it is YES, go to STEPs 6 & 7, otherwise go    //
        //                directly to STEP 8.                                      //
        //           6. Copy the secondary level set function to the primary levle //
        //                set function through comparision between these two and   //
        //                taking the minimal value for the primarary lelel set     //
        //           7. Use the Fast Marching Method to stretch this updated       //
        //                primary level set function                               //
        //           8. Update primary level set function by (1) extending boundary//
        //                points' velocities to narrow band points, (2) computing  //
        //                gradient of primary level set function, (3) using the    //
        //                Upwind Difference method to update the primary level set //
        //           9. Check convergence criterion. If YES, terminating the       //
        //                program, otherwise repeating STEPs 1-9 and other STEPs   //
        //                relating to primary level set function.                  //
        //                                                                         //
        /////////////////////////////////////////////////////////////////////////////
        std::vector<double> signedDistance_temp(lsm_mesh.nNodes);
        if (n_iterations > 5) {

        	if (isHole) {

	            //
	            // Step 1.
	            //

	            h_count = hole_map(lsm_mesh, level_set, h, lBand, h_index, h_elem);

	            //
	            // Step 3. Compute node sensitivities using least squares method
	            //

	            // 3.1 Extrapolate nodal sensitvities from sensitivites for gauss points
	            for (int inode = 0; inode < lsm_mesh.nNodes; inode++ ) {
	                vector<double> nPoint (2, 0);
	                nPoint[0] = lsm_mesh.nodes[inode].coord.x;
	                nPoint[1] = lsm_mesh.nodes[inode].coord.y;

	                // Interpolate Node Point sensitivities by least squares.
	                // sens.ComputeBoundarySensitivities(radius, nPoint) ;
	                sens.ComputeBoundarySensitivities(nPoint) ;

	                // Assign sensitivities.
	                h_Nsens_Temp[inode].sensitivities.resize(2);
	                fill(h_Nsens_Temp[inode].sensitivities.begin(), h_Nsens_Temp[inode].sensitivities.end(), 0.0);

	                h_Nsens_Temp[inode].sensitivities[0] = -sens.boundary_sensitivities[inode] ;
	                h_Nsens_Temp[inode].sensitivities[1] = -1 ;
	            }
	            // clearing sens.boundarysens vector.
	            sens.boundary_sensitivities.clear() ;

	            // 3.2 Calculate element sensitivities
	            for (int iel = 0; iel < lsm_mesh.nElements; iel++ ) {

	                h_Esens[iel].sensitivities.resize(2);
	                fill(h_Esens[iel].sensitivities.begin(), h_Esens[iel].sensitivities.end(), 0.0);
	                for (int ind = 0; ind < 4; ind++) {
	                    int inode;
	                    inode = lsm_mesh.elements[iel].nodes[ind];
	                    h_Esens[iel].sensitivities[0] += 0.25 * h_Nsens_Temp[inode].sensitivities[0];
	                    h_Esens[iel].sensitivities[1] += 0.25 * h_Nsens_Temp[inode].sensitivities[1];
	                }

	            }
	            // 3.3 Update nodal sensitivities
	            // clear nodal sensitivity value
	            for (int inode = 0; inode < lsm_mesh.nNodes; inode++ ) {
	                h_Nsens[inode].sensitivities.resize(2);
	                fill(h_Nsens[inode].sensitivities.begin(), h_Nsens[inode].sensitivities.end(), 0.0);
	            }
	            // re-assign nodal sensitivity value based on calculatd element sensitivities
	            for (int iel = 0; iel < lsm_mesh.nElements; iel++ ) {
	                for (int ind = 0; ind < 4; ind++) {
	                    int inode;
	                    inode = lsm_mesh.elements[iel].nodes[ind];
	                    h_Nsens[inode].sensitivities[0] += 0.25 * h_Esens[iel].sensitivities[0];
	                    h_Nsens[inode].sensitivities[1] += 0.25 * h_Esens[iel].sensitivities[1];
	                }
	            }

	            if (h_flag) {
	                //
	                // Step 2. Initialise the secondary level set function after inserting new holes
	                //
	                h_lsf.resize(lsm_mesh.nNodes); fill(h_lsf.begin(), h_lsf.end(), h_bar);//}
	                get_h_lsf( lsm_mesh.nNodes, h_index, h_Nsens, lambdas, h_lsf );
	                //initialise_hole_lsf(lsm_mesh, h_count, holeCFL, level_set, moveLimit, h_index, h_elem, h_Nsens, h_lsf, lambdas);
	                h_flag = false;
	            } else {
	                //
	                // Step 4. Update existing h_lsf when new hole has not yet been activated.
	                //
	                get_h_lsf( lsm_mesh.nNodes, h_index, h_Nsens, lambdas, h_lsf );
	                //
	                // Step 5. Check whether new holes should be inserted
	                //
	                int inserted_hole_nodes = 0;
	                for (int inode = 0; inode < lsm_mesh.nNodes; inode++ ) {
	                    if ( (h_index[inode] ==1) && (h_lsf[inode] < 0) ) {
	                        if (!h_flag) {
	                            cout << "\n\n--------------------------------------------\n";
	                            cout << "  Hole will be inserted at ";
	                        }
	                        h_flag = true;
	                        inserted_hole_nodes ++;
	                        //cout<< "\n" << inode << "\t" << h_lsf[inode] << "\t" << level_set.signedDistance[inode];
	                     }
	                }
	                // counting hole_area secondary
	                double area_h_lsf, area_lsf;
	                area_h_lsf = LSM::Boundary_hole(level_set, &h_lsf).computeAreaFractions();
	                area_lsf   = LSM::Boundary_hole(level_set, &level_set.signedDistance).computeAreaFractions();

	                printf("\nThe area fraction corresponding to lsf and h_lsf are: %8.2f %8.2f %8.2f\n", boundary.area, area_lsf, area_h_lsf) ;

	                //cout << "\n\nNumber of nodes to be set as new hole nodes: " << inserted_hole_nodes << endl;

	                //
	                // Step 6. Copy values of secondary level set function to primaray level set function
	                //
	                io.saveLevelSetVTK(9000, level_set) ;

	                if (h_flag) {

	                    //
	                    // 6.1 Check whether new holes' area exceeds a certain threshold
	                    //

	                    double hole_areafraction = (mesh_area- area_h_lsf)/area_lsf;
	                    // Find minimum of h_lsf
	                    double temp_min_h_lsf = 1.0;
	                    if (hole_areafraction > newHoleAreaLimit) {
	                        for (int inode = 0; inode < lsm_mesh.nNodes; inode++) {
	                            temp_min_h_lsf = (temp_min_h_lsf < h_lsf[inode]) ? temp_min_h_lsf : h_lsf[inode];
	                        }
	                    }
	                    while (hole_areafraction > newHoleAreaLimit) {

	                        cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	                        cout << "     Too much material is removed. \n";
	                        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

	                        // move up h_lsf untill the limit of the area of new holes to be inserted meeting the requirement
	                        for (int inode = 0; inode < lsm_mesh.nNodes; inode++) {
	                            h_lsf[inode] = h_lsf[inode] - 0.005*temp_min_h_lsf;
	                        }
	                        area_h_lsf = LSM::Boundary_hole(level_set, &h_lsf).computeAreaFractions();
	                        hole_areafraction = (mesh_area- area_h_lsf)/area_lsf;
	                    }

	                    //
	                    // 6.2 Update primary level set function to insert new holes
	                    //
	                    signedDistance_temp.clear(); double min_h_lsf = 1.0, min_lsf = 1.0;
	                    for (int inode = 0; inode < lsm_mesh.nNodes; inode++) {
	                        // signedDistance_temp[inode] = level_set.signedDistance[inode];
	                        if ( (h_index[inode] ==1) && (hole_areafraction>1e-3)) {
	                            //if (h_lsf[inode] < level_set.signedDistance[inode]) {
	                             if ((h_lsf[inode]<=h_bar) && (h_lsf[inode] < level_set.signedDistance[inode])) {
	                                level_set.signedDistance[inode] = h_lsf[inode];
	                            }
	                        }
	                        signedDistance_temp[inode] = level_set.signedDistance[inode];
	                        min_h_lsf = (min_h_lsf < h_lsf[inode]) ? min_h_lsf : h_lsf[inode];
	                        min_lsf = (min_lsf < level_set.signedDistance[inode]) ? min_lsf : level_set.signedDistance[inode];
	                    }
	                    cout <<"\n\nMininal primary and secondary LSF: " << min_h_lsf << "\t" << min_lsf << endl;
	                }
	                //
	                // Step 7. Use fast marching method to re-initialise signed distance function
	                //
	                if (h_flag) {

	                    // for (int inode = 0; inode < lsm_mesh.nNodes; inode++ ) {
	                    //     level_set.signedDistance[inode] = h_lsf[inode];
	                    // }
	                    // cout << "\nThe area fraction corresponding to h_lsf is [new] : " << LSM::Boundary_hole(level_set,&h_lsf).computeAreaFractions() << endl;
	                    // io.savelevel_setVTK(9001, level_set) ;

	                    // for (int inode = 0; inode < lsm_mesh.nNodes; inode++ ) {
	                    //     level_set.signedDistance[inode] = signedDistance_temp[inode];
	                    // }
	                    // io.savelevel_setVTK(9002, level_set) ;

	                    // // boundary_hole2.computeAreaFractions() ;
	                    // cout << "\nThe area fraction corresponding to lsf before stretching is: " << LSM::Boundary_hole(level_set,&level_set.signedDistance).computeAreaFractions() << endl;

	                    M2DO_LSM::FastMarchingMethod fmm(lsm_mesh, false);
	                    fmm.march(level_set.signedDistance);

	                    // io.savelevel_setVTK(9003, level_set) ;

	                    // cout << "\nhe area fraction corresponding to lsf after stretching is: " << LSM::Boundary_hole(level_set,&level_set.signedDistance).computeAreaFractions() << endl;

	                }

	            }

	            if (h_flag) { cout << "\n--------------------------------------------\n\n"; }
	        }

        }


		////////////////////////////////////////////////////////////////////////////

		// Extend boundary point velocities to all narrow band nodes
		level_set.computeVelocities (boundary.points, time_step, 0, rng) ;

		// Compute gradient of the signed distance function within the narrow band
		level_set.computeGradients () ;

		// Update the level set function
		bool is_reinitialised = level_set.update (time_step) ;

		// Reinitialise the signed distance function, if necessary
		if (!is_reinitialised) {
			// Reinitialise at least every 20 iterations
			if (n_reinit == 20) {
				level_set.reinitialise () ;
				n_reinit = 0 ;
			}

		} else n_reinit = 0 ;

		// Increment the number of steps since reinitialisation
		n_reinit++ ;

		// Increment the time
		time += time_step ;

		// Calculate current area fraction
		double area = boundary.area / mesh_area ;

		// Record the time and area
		times.push_back (time) ;
		areas.push_back (area) ;

		// Converence criterion [Dunning_11_FINAL]:
		// find the max relative distance over the past five iterations:
		objective_values.push_back (sens.objective) ;
		double objective_value_k, objective_value_m ;

		if (n_iterations > 5) {

			objective_value_k = sens.objective ;
			relative_difference = 0.0 ;

			for (int i = 1 ; i <= 5 ; i++) {
				objective_value_m = objective_values[n_iterations - i - 1] ;
				relative_difference = max(relative_difference, abs((objective_value_k - objective_value_m)/objective_value_k)) ;
			}

		}

		// Print statistics
		printf ("%8.1f %12.4f %10.4f\n", double(n_iterations), sens.objective, area) ;

		// Print statistics to .txt file
		history_file.open ("results/history/history.txt", ios_base::app) ;
		history_file << n_iterations << "\t" << sens.objective << "\t" << area << "\n" ;
		history_file.close () ;

		// Write level set and area fractions to .vtk file
		io.saveLevelSetVTK (n_iterations, level_set, false, false, "results/level_set") ;
		io.saveAreaFractionsVTK (n_iterations, lsm_mesh, "results/area_fractions") ;

		// Write level set, area fractions, and boundary segments to .txt file:
		io.saveBoundarySegmentsTXT (n_iterations, boundary, "results/boundary_segments") ;

		// Check if convergence has been met:
		if ((relative_difference < max_diff) & (area < 1.001 * max_area)) break;

	}

	// END OF LEVEL SET TOPOLOGY OPTIMIZATION LOOP


	cout << "\nProgram complete.\n" << flush ;
	return 0 ;

}

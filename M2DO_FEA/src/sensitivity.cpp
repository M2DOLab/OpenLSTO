#include "stationary_study.h"
#include "sensitivity.h"

using namespace M2DO_FEA ;

// Constructor
SensitivityAnalysis :: SensitivityAnalysis (StationaryStudy & study) : study (study) {

    order    = study.mesh.spacedim ;
    spacedim = study.mesh.spacedim ;

	// Declaring variables
	int number_of_elements = study.mesh.solid_elements.size();   // Total number of elements
	int number_of_gauss_points = pow(order, spacedim); // Total number of gauss points

	// Resizing sensitivities in function of number_of_elements and number_of_gauss_points
	sensitivities.resize(number_of_elements);
	for (int i = 0; i < number_of_elements; i++)
	{
		sensitivities[i].sensitivity_at_gauss_point.resize(number_of_gauss_points);
		sensitivities[i].sensitivity_component1_at_gauss_point.resize(number_of_gauss_points);
		sensitivities[i].sensitivity_component2_at_gauss_point.resize(number_of_gauss_points);
		sensitivities[i].sensitivity_component3_at_gauss_point.resize(number_of_gauss_points);
		sensitivities[i].coordinate.resize(number_of_gauss_points);
		for (int j = 0; j < number_of_gauss_points; j++)
		{
			sensitivities[i].coordinate[j].resize(spacedim);
		}
	}

    // Computing sensitivity coordinates.
    ComputeSensitivitiesCoordinates (false) ;
    // Compute element centroids.
    // study.mesh.ComputeCentroids () ;
}


// Functionality to compute compliance sensitivities for use in LSM

void SensitivityAnalysis :: ComputeSensitivitiesCoordinates (bool time_it) {

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing sensitivity coordinates ... " << flush ;

    }

	// Declaring variables
	int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
	int number_of_gauss_points = pow(order,spacedim); // Total number of gauss points
	vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
	// Shape functions and quadrature objects
	LinearShapeFunction linear_shape_function (spacedim, spacedim) ;
    GaussianQuadrature  quadrature (spacedim, order) ;

	// For each element i
    for (int i = 0; i < number_of_elements; i++)
    {
    	// For each gauss point
    	for (int j = 0; j < number_of_gauss_points; j++)
    	{
    		// Selecting gauss points (in order to compute their global coordinates)
        	for (int k = 0 ; k < spacedim; k++)
        	{
				eta[k]  = quadrature.eta[eta_count[k]];
			}

    		// Storing gauss point coordinates
    		for (int k = 0; k < spacedim; k++)
    		{
    			sensitivities[i].coordinate[j][k] = 0.0;
				for (int l = 0; l < pow(2, spacedim); l++)
				{
					sensitivities[i].coordinate[j][k] += linear_shape_function.GetShapeFunctionValues(l, eta)*study.mesh.nodes[study.mesh.solid_elements[i].node_ids[l]].coordinates[k];
				}
    		}

    		// Update eta counter (to select next group of Gauss points)
			eta_count = quadrature.UpdateEtaCounter(eta_count);
    	}
    }

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }

}

// COMPLIANCE PROBLEM //

// Functionality to compute compliance sensitivities for use in LSM

void SensitivityAnalysis :: ComputeComplianceSensitivities (bool time_it) {

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing compliance sensitivities ... " << flush ;

    }

	// DECLARING VARIABLES AND RESIZING VECTORS

	// Scalars and vectors
	int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
	int number_of_gauss_points = pow(order,spacedim); // Total number of gauss points
	vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
	VectorXd element_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of element displacements
	vector<int> dof; // Vector with dofs

	// Stress*strain matrices.
    vector<MatrixXd> B;
    B.resize(number_of_gauss_points);
	MatrixXd Bu = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd C = study.mesh.solid_materials[0].C ;
    MatrixXd stress_strain;

	// Quadrature object.
    GaussianQuadrature  quadrature (spacedim, order) ;

	// FUNCTION BODY

    // Computing strain-displacement matrices.
    for (int j = 0; j < number_of_gauss_points; j++)
    {
        // Selecting Gauss points.
        for (int k = 0 ; k < spacedim; k++)
        {
            eta[k]  = quadrature.eta[eta_count[k]];
        }

        // Strain-displacement matrix at Gauss point.
        B[j] = study.mesh.solid_elements[0].B(eta);

        // Update eta counter (to select next group of Gauss points).
        eta_count = quadrature.UpdateEtaCounter(eta_count);
    }

	// For each element i
    for (int i = 0; i < number_of_elements; i++)
    {
        // If the element is too soft (very small area fraction)
        if (study.mesh.solid_elements[i].area_fraction <= 0.1)
        {
        	// For each gauss point
        	for (int j = 0; j < number_of_gauss_points; j++)
        	{
        		// Sensitivity is not computed and set as zero
        		sensitivities[i].sensitivity_at_gauss_point[j] = 0.0;
        	}
        }
        // If the element has significant area fraction
        else
        {
        	// For each Gauss point
        	for (int j = 0; j < number_of_gauss_points; j++)
        	{
				// Element dofs
                dof = study.mesh.solid_elements[i].dof ;

				// Selecting element displacements.
				for (int k = 0 ; k < dof.size() ; k++)
				{
					element_displacements(k) = study.u(dof[k]) ;
				}

				// Strain.
                Bu = B[j]*element_displacements;

				// Sensitivities (stress*strain).
                stress_strain = Bu.transpose()*C*Bu;
                sensitivities[i].sensitivity_at_gauss_point[j] = -stress_strain(0,0)*(study.mesh.solid_elements[i].area_fraction);

        	}
        }
    }

    // Objective function (compliance).
    objective = study.f.dot(study.u);

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }

}


// STRESS PROBLEM //

// P-norm stress shape sensitivities.
void SensitivityAnalysis :: ComputeStressSensitivities (bool time_it, double pnorm) {

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing stress sensitivities ... " << flush ;

    }

    // DECLARING VARIABLES AND RESIZING VECTORS

    // Scalars and vectors
    int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
    int number_of_gauss_points = pow(2,spacedim); // Total number of gauss points
    vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
    VectorXd element_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of element displacements
    VectorXd adjoint_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of adjoint displacements
    vector<int> dof; // Vector with dofs

    // Stress*strain matrices.
    vector<MatrixXd> B;
    B.resize(number_of_gauss_points);
    MatrixXd Bu = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd Bu_adj = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd CBu = MatrixXd::Zero(spacedim*spacedim, 1);
    MatrixXd CB = Eigen::MatrixXd::Zero(pow(spacedim,spacedim), pow(2,spacedim)*spacedim); // and C*B
    MatrixXd stress_strain_adj;
    MatrixXd Tvm2, lambda_i;
    double Tvm;

    // Quadrature object.
    GaussianQuadrature  quadrature (spacedim, order) ;

    // Resizing strain vector.
    strains.resize(number_of_elements);
    for (int i = 0; i < number_of_elements; i++)
    {
        strains[i].von_mises.resize(number_of_gauss_points);
        for (int j = 0; j < number_of_gauss_points; j++)
        {
            strains[i].von_mises[j] = 0.0;
        }
    }



    // Initializing adjoint vector with zeros (then it will be assembled for each element pseudo-load).
    int n_dof = study.dirichlet_boundary_conditions.mesh_n_dof;
    study.f_i = VectorXd::Zero(n_dof);
    // Number of dofs for reduced matrix (after Dirichlet BC's).
    int n_dof_reduced = study.dirichlet_boundary_conditions.mesh_n_dof - study.dirichlet_boundary_conditions.dof.size();
    study.f_i_reduced = VectorXd::Zero(n_dof_reduced); // Initializing adjoint pseudo-load vector.

    // FUNCTION BODY

    // Computing strain-displacement matrices.
    for (int j = 0; j < number_of_gauss_points; j++)
    {
        // Selecting Gauss points.
        for (int k = 0 ; k < spacedim; k++)
        {
            eta[k]  = quadrature.eta[eta_count[k]];
        }

        // Strain-displacement matrix at Gauss point.
        B[j] = study.mesh.solid_elements[0].B(eta);

        // Update eta counter (to select next group of Gauss points).
        eta_count = quadrature.UpdateEtaCounter(eta_count);
    }

    // Integral of B.
    MatrixXd B_int = study.mesh.solid_elements[0].B_int();

    // Elasticity tensor.
    MatrixXd C = study.mesh.solid_materials[0].C ;

    // Voigt matrix.
    MatrixXd Voigt = study.mesh.solid_materials[0].V;

    // Element volume.
    // double volume = study.mesh.solid_elements[0].V();

    // Initializing objective.
    objective = 0.0;

    // Initializing maximum stress.
    von_mises_max = 0.0;

    // For each element (compute von Mises stress and pseudo-load).
    for (int i = 0; i < number_of_elements; i++)
    {
        // Check if the element is out of the domain or not.
        if (study.mesh.solid_elements[i].area_fraction <= 0.1) // Out.
        {
            // Do nothing.
        }
        else // In.
        {
            // Element dofs
            dof = study.mesh.dof(study.mesh.solid_elements[i].node_ids);

            // Selecting element displacements.
            for (int j = 0 ; j < dof.size() ; j++)
            {
                element_displacements(j) = study.u(dof[j]);
            }

            // C*B*u and C*B at point eta.
            CBu = C*B_int*element_displacements;
            CB = C*B_int;

            // von Mises stress.
            Tvm2 = CBu.transpose()*Voigt*CBu;
            Tvm = sqrt(Tvm2(0,0));

            if (!sensitivities[i].isExcluded)
            {
                // Sum stress to the p-norm for non-excluded elements.
                objective += pow(study.mesh.solid_elements[i].area_fraction*Tvm,pnorm);

                // Storing maximum von Mises stress in the domain.
                if (study.mesh.solid_elements[i].area_fraction*Tvm > von_mises_max)
                { von_mises_max = study.mesh.solid_elements[i].area_fraction*Tvm; }
            }


            // Adjoint force vector.
            lambda_i = -pnorm*pow(Tvm,pnorm-2)*CB.transpose()*Voigt*CBu;

            // Asemble adjoint force vector.
            study.AssembleF_i(lambda_i, dof, false);
        }
    }

    // Objective (p-norm function).
    objective = pow(objective,1/pnorm);

    // Solve adjoint equation.
    // study.SolveWithHSLMA57_f_i (false, false, false) ;
    study.SolveWithCG_f_i();

    // For each element i
    for (int i = 0; i < number_of_elements; i++)
    {
        // Initializing average measures.
        sensitivities[i].sensitivity_average = 0.0; //Sensitivities
        strains[i].von_mises_average = 0.0; // Stresses.

        // If the element is too soft (very small area fraction)
        if (study.mesh.solid_elements[i].area_fraction <= 0.1)
        {
            // For each gauss point
            for (int j = 0; j < number_of_gauss_points; j++)
            {
                // Sensitivity is not computed and set as zero
                sensitivities[i].sensitivity_at_gauss_point[j] = 0.0;
		sensitivities[i].sensitivity_component1_at_gauss_point[j] = 0.0;
		sensitivities[i].sensitivity_component2_at_gauss_point[j] = 0.0;
            }
        }
        // If the element has significant area fraction
        else
        {
            // Element dofs
            dof = study.mesh.dof(study.mesh.solid_elements[i].node_ids) ;

            // Selecting element displacements
            for (int j = 0 ; j < dof.size(); j++)
            {
                element_displacements(j) = study.u(dof[j]) ;
                adjoint_displacements(j) = study.u_i(dof[j]);
            }

            // For each Gauss point
            for (int j = 0; j < number_of_gauss_points; j++)
            {
                // C*B*u at point eta.
                CBu = C*B[j]*element_displacements;
                Bu_adj = B[j]*adjoint_displacements;

                // Stress(mechanical)*strain(adjoint)
                MatrixXd stress_strain_adj = CBu.transpose()*Bu_adj;

                // von Mises stress.
                Tvm2 = CBu.transpose()*Voigt*CBu;
                strains[i].von_mises[j] = sqrt(Tvm2(0,0));

                // Storing stress average per element.
                strains[i].von_mises_average += strains[i].von_mises[j]/number_of_gauss_points;

                // Sensitivities (stress*strain).
                sensitivities[i].sensitivity_at_gauss_point[j] = pow(objective,1-pnorm)*(pow(strains[i].von_mises[j],pnorm) + stress_strain_adj(0,0))/pnorm;

		// Sensitivities (stress*strain).
		// Components of sensitivity (stress*strain)
		sensitivities[i].sensitivity_component1_at_gauss_point[j] = strains[i].von_mises[j] * study.mesh.solid_elements[i].area_fraction;
		sensitivities[i].sensitivity_component2_at_gauss_point[j] = stress_strain_adj(0,0) * study.mesh.solid_elements[i].area_fraction;

		// Storing sensitivity average per element
		sensitivities[i].sensitivity_average += sensitivities[i].sensitivity_at_gauss_point[j]/number_of_gauss_points;
            }
        }
    }

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }
}

void SensitivityAnalysis :: ComputeStressSensitivities3D (bool time_it, double pnorm ) {

  //Same as 2D stress sensitivities, except for 3d
    double qnorm = pnorm;
    double area_min = 0.0;

    auto t_start = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "\nComputing stress sensitivities ... " << flush ;

    }

    // DECLARING VARIABLES AND RESIZING VECTORS

    // Scalars and vectors
    int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
    int number_of_gauss_points = pow(order,spacedim); // Total number of gauss points
    vector<double> eta(spacedim,0), eta_count(spacedim,0); // Vector of gauss points
    VectorXd element_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of element displacements
    VectorXd adjoint_displacements = VectorXd::Zero(pow(2,spacedim)*spacedim); // Vector of adjoint displacements
    vector<int> dof; // Vector with dofs

    // Stress*strain matrices.
    vector<MatrixXd> B;
    B.resize(number_of_gauss_points);
    MatrixXd Bu = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd Bu_adj = MatrixXd::Zero(pow(spacedim,spacedim), 1);
    MatrixXd CBu = MatrixXd::Zero(spacedim*spacedim, 1);
    MatrixXd CB = Eigen::MatrixXd::Zero(pow(spacedim,spacedim), pow(2,spacedim)*spacedim); // and C*B
    MatrixXd stress_strain_adj;
    MatrixXd Tvm2, lambda_i;
    double Tvm;

    // Integral of B.
    MatrixXd B_int = study.mesh.solid_elements[0].B_int();

    // Resizing strain vector.
    strains.resize(number_of_elements);

    // Quadrature object.
    GaussianQuadrature  quadrature (spacedim, order) ;

    // Initializing adjoint vector with zeros (then it will be assembled for each element pseudo-load).
    int n_dof = study.dirichlet_boundary_conditions.mesh_n_dof;
    study.f_i = VectorXd::Zero(n_dof);
    // Number of dofs for reduced matrix (after Dirichlet BC's).
    int n_dof_reduced = study.dirichlet_boundary_conditions.mesh_n_dof - study.dirichlet_boundary_conditions.dof.size();
    study.f_i_reduced = VectorXd::Zero(n_dof_reduced); // Initializing adjoint pseudo-load vector.


    // FUNCTION BODY

    // Computing strain-displacement matrices.
    for (int j = 0; j < number_of_gauss_points; j++)
    {
        // Selecting Gauss points.
        for (int k = 0 ; k < spacedim; k++)
        {
            eta[k]  = quadrature.eta[eta_count[k]];
        }

        // Strain-displacement matrix at Gauss point.
        B[j] = study.mesh.solid_elements[0].B(eta);

        // Update eta counter (to select next group of Gauss points).
        eta_count = quadrature.UpdateEtaCounter(eta_count);
    }


    // Elasticity tensor.
    MatrixXd C = study.mesh.solid_materials[0].C ;

    // Voigt matrix.
    MatrixXd Voigt = study.mesh.solid_materials[0].V;

    // Element volume.
    // double volume = study.mesh.solid_elements[0].V();

    // Initializing objective.
    objective = 0.0;

    // Initializing maximum stress.
    von_mises_max = 0.0;

    // For each element (compute von Mises stress and pseudo-load).
    for (int i = 0; i < number_of_elements; i++)
    {
        // Check if the element is out of the domain or not.
        if (study.mesh.solid_elements[i].area_fraction < area_min) // Out.
        {
            // Do nothing.
        }
        else // In.
        {
            // Element dofs
            dof = study.mesh.dof(study.mesh.solid_elements[i].node_ids);

            // Selecting element displacements.
            for (int j = 0 ; j < dof.size() ; j++)
            {
                element_displacements(j) = study.u(dof[j]);
            }

            MatrixXd Tvmmat_avg = (C * (B_int *element_displacements ) ).transpose() * Voigt * (  C * (B_int *element_displacements )  );



            // for each gauss point
            for (int j = 0; j < number_of_gauss_points; j++)
            {
              // Define Tvm
              MatrixXd Tvmmat = (C * (B[j] *element_displacements ) ).transpose() * Voigt * (  C * (B[j] *element_displacements )  );
              Tvm = sqrt(Tvmmat(0,0));

              // Adjoint force vector.
              lambda_i = pow(study.mesh.solid_elements[i].area_fraction  , qnorm )*pnorm*pow(Tvm,pnorm-2)* ( (C*B[j]).transpose() * Voigt * C) *(B[j] *element_displacements );

              // Asemble adjoint force vector.
              study.AssembleF_i(lambda_i, dof, false);

              objective += pow(study.mesh.solid_elements[i].area_fraction  , qnorm )*pow(Tvm , pnorm )/number_of_gauss_points ;

              if ( Tvm*study.mesh.solid_elements[i].area_fraction > von_mises_max)
              { von_mises_max = Tvm*study.mesh.solid_elements[i].area_fraction ; }


            }

        }
    }



    // Objective (p-norm function).
    objective = pow(objective,1.0/pnorm);

    // Solve adjoint equation.
    study.SolveWithCG_f_i () ;

    // For each element i
    for (int i = 0; i < number_of_elements; i++)
    {

      strains[i].von_mises_average = 0.0; // Stresses.

        // If the element is too soft (very small area fraction)
        if (study.mesh.solid_elements[i].area_fraction < area_min)
        {
            // For each gauss point
            for (int j = 0; j < number_of_gauss_points; j++)
            {
                // Sensitivity is not computed and set as zero
                sensitivities[i].sensitivity_at_gauss_point[j] = 0.0;

                // components of the sensitivities
                sensitivities[i].sensitivity_component1_at_gauss_point[j] = 0.0;
                sensitivities[i].sensitivity_component2_at_gauss_point[j] = 0.0 ;
                sensitivities[i].sensitivity_component3_at_gauss_point[j] = 0.0 ;
            }
        }
        // If the element has significant area fraction
        else
        {
            // Element dofs
            dof = study.mesh.dof(study.mesh.solid_elements[i].node_ids) ;

            // Selecting element displacements
            for (int j = 0 ; j < dof.size(); j++)
            {
                element_displacements(j) = study.u(dof[j]) ;
                adjoint_displacements(j) = study.u_i(dof[j]);
            }

            // For each Gauss point
            for (int j = 0; j < number_of_gauss_points; j++)
            {
                // C*B*u at point eta.
                CBu = C*B[j]*element_displacements;
                Bu_adj = B[j]*adjoint_displacements;

                // Stress(mechanical)*strain(adjoint)
                MatrixXd stress_strain_adj = CBu.transpose()*Bu_adj;

                // Define Tvm
                MatrixXd Tvmmat = (C * (B[j] *element_displacements) ).transpose() * Voigt * (  C * (B[j] *element_displacements )  );
                Tvm = sqrt(Tvmmat(0,0));

                // Sensitivities (stress*strain).
                MatrixXd stress_strain_mat_me = Bu_adj.transpose()*CBu ;

                // sensitivities at the gauss point
                sensitivities[i].sensitivity_at_gauss_point[j] = - stress_strain_mat_me(0,0) +
                pow(Tvm , pnorm)*pow(study.mesh.solid_elements[i].area_fraction  , qnorm - 1)*qnorm;

                // components of the sensitivities
                sensitivities[i].sensitivity_component1_at_gauss_point[j] = Tvm * pow(study.mesh.solid_elements[i].area_fraction  , 1.0*(qnorm - 1)/(1.0*pnorm)) * pow(qnorm , 1.0/(1.0*pnorm) );
                sensitivities[i].sensitivity_component2_at_gauss_point[j] = stress_strain_mat_me(0,0) ;

                strains[i].von_mises_average += study.mesh.solid_elements[i].area_fraction*Tvm/number_of_gauss_points; // Stresses.

            }
        }
    }

    auto t_end = chrono::high_resolution_clock::now() ;

    if (time_it) {

        cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

    }
}



// LEAST SQUARES FUNCTIONALITIES

// Functionality to compute sensitivities at the level-set boundaries.

void SensitivityAnalysis :: ComputeBoundarySensitivities (vector<double> boundary_point, double radius, int indicator, double p_norm) {

    // DECLARING VARIABLES

    // Scalars and vectors
    int number_of_elements = study.mesh.solid_elements.size(); // Total number of elements
    int number_of_gauss_points = pow(order, spacedim); // Total number of gauss points
    double squared_radius = radius*radius; // Squared radius.
    double distance; // Squared distance from the Gauss point to the element i.

    // Least square information object.
    LeastSquares leastsq;

    // FUNCTION BODY

    // Finding Gauss points inside the boundary and the subdomain defined by radius
    for (int i = 0; i < number_of_elements; i++)
    {
        // Looking at nearby elements
        int aux_out = 0;
        for (int j = 0; j < spacedim; j++)
        {
            // If the element is too far from boundary point -> aux_out = 1.
            if ((study.mesh.solid_elements[i].centroid[j] > (boundary_point[j]+1.5*radius)) || (study.mesh.solid_elements[i].centroid[j] < (boundary_point[j]-1.5*radius)))
            {
                aux_out = 1;
            }
        }
        // If element is outside the subdomain defined by 1.5radius
        if (aux_out == 1)
        {
            // Element is too far from boundary point
        }
        else // Else element is close enough to be considered
        {
            // Check if element is inside the boundary
            if (study.mesh.solid_elements[i].area_fraction >= 0.0)
            {
                // For each Gauss point of the ith element
                for (int j = 0; j < number_of_gauss_points; j++)
                {
                    // Compute squared distance between current Gauss point and boundary point
                    distance = 0.0;
                    for (int k = 0; k < spacedim; k++)
                    {
                        // Squared distance in each space dimension
                        distance += pow(boundary_point[k]-sensitivities[i].coordinate[j][k],2);
                    }

                    // If squared distance is less than squared radius, save information
                    if (distance < squared_radius)
                    {
                        // Grouping information in object
                        leastsq.distance_from_gauss_point = sqrt(distance);      // Save distance
                        leastsq.area_fraction_at_gauss_point = study.mesh.solid_elements[i].area_fraction;  // Save area fraction
                        leastsq.element_number = i;     // Save element number
                        leastsq.gauss_point_number = j;  // Save Gauss point number
                        leastsq.coordinate = sensitivities[i].coordinate[j];   // Coordinates
                        // Storing information in class
                        least_squares.push_back(leastsq);
                    }
                }
            }
        }
    }

    // If not enough Gauss points inside subdomain, we are on an island.
    if (least_squares.size() < 10)
    {
        // Sensitivities at the boundaries for islands are zero.
        boundary_sensitivities.push_back(0.0);
        return;
    }

    // SOLVING LEAST SQUARES
    double B;

    // Solve least squares and obtain sensitivity at the boundary
    switch (indicator) {

    case 0 : { // compliance (default) case
      B = SolveLeastSquares(least_squares, boundary_point);
      break;
    }

    case 1 : { // stress case
      double B1 = SolveLeastSquares(least_squares, boundary_point, 1);
      double B2 = SolveLeastSquares(least_squares, boundary_point, 2);
      B = pow(objective, 1 - p_norm) * (pow(B1, p_norm) + B2) / p_norm;
      break;
    }

    }

    // Store sensitivity
    boundary_sensitivities.push_back(B);

    // Clear leastsquares vector.
    least_squares.clear();

}

// Functionality to solve least squares problem for 2D and 3D cases.

double SensitivityAnalysis :: SolveLeastSquares(vector<LeastSquares> least_squares, vector<double> boundary_point, int indicator) {

    // Number of Gauss points inside subdomain defined by radius.
    int number_of_gauss_points = least_squares.size();

    // Size of elements in least squares' basis function.
    int basis[2] = {6,10};
    // Selecting size according to dimensionality of the problem.
    int n = basis[spacedim-2];

    //initialize sensitivity
    double sens;

    std::vector<std::vector<double>> A_dash(number_of_gauss_points);
    for(int i = 0; i< number_of_gauss_points; i++) A_dash[i].resize(n,0.0);

    std::vector<double> B_dash(number_of_gauss_points,0.0);

    // Building least squares problem.
    if (spacedim == 2) // For 2D interpolation.
    {
        for (int i = 0; i < number_of_gauss_points; i++)
        {
            // Weight function by inverse distance.
            double lsweight = least_squares[i].area_fraction_at_gauss_point/least_squares[i].distance_from_gauss_point;

            // Relative x and y coordinates.
            double xb = least_squares[i].coordinate[0] - boundary_point[0];
            double yb = least_squares[i].coordinate[1] - boundary_point[1];

            // Storing weighted distances information.

            A_dash[i][0] = lsweight;
            A_dash[i][1] = xb * lsweight;
            A_dash[i][2] = yb * lsweight;
            A_dash[i][3] = xb * yb * lsweight;
            A_dash[i][4] = xb * xb * lsweight;
            A_dash[i][5] = yb * yb * lsweight;



            // Sensitivity at the current point.
	    switch (indicator) {

	    case 0 : //default case, solve for sensitivity
	      sens = sensitivities[least_squares[i].element_number].sensitivity_at_gauss_point[least_squares[i].gauss_point_number];
	      break;

	    case 1 : // solve for first component of sensitivity
	      sens = sensitivities[least_squares[i].element_number].sensitivity_component1_at_gauss_point[least_squares[i].gauss_point_number];
	      break;

	    case 2 : // solve for second component of sensitivity
	      sens = sensitivities[least_squares[i].element_number].sensitivity_component2_at_gauss_point[least_squares[i].gauss_point_number];
	      break;

	    case 3 : // solve for third component of sensitivity
	      sens = sensitivities[least_squares[i].element_number].sensitivity_component3_at_gauss_point[least_squares[i].gauss_point_number];
	      break;

	    }

            // Storing weighted sensitivity.

            B_dash[i] = sens*lsweight;
        }
    }
    else if (spacedim == 3) // For 3D interpolation.
    {
        for (int i = 0; i < number_of_gauss_points; i++)
        {
            // Weight function by inverse distance.
            double lsweight = least_squares[i].area_fraction_at_gauss_point/least_squares[i].distance_from_gauss_point;

            // Relative x, y and z coordinates.
            double xb = least_squares[i].coordinate[0] - boundary_point[0];
            double yb = least_squares[i].coordinate[1] - boundary_point[1];
            double zb = least_squares[i].coordinate[2] - boundary_point[2];

            // Storing weighted distances information.
            A_dash[i][0] = lsweight;
            A_dash[i][1] = xb * lsweight;
            A_dash[i][2] = yb * lsweight;
            A_dash[i][3] = zb * lsweight;

            A_dash[i][4] = xb * yb * lsweight;
            A_dash[i][5] = xb * zb * lsweight;
            A_dash[i][6] = yb * zb * lsweight;

            A_dash[i][7] = xb * xb * lsweight;
            A_dash[i][8] = yb * yb * lsweight;
            A_dash[i][9] = zb * zb * lsweight;

            switch (indicator) {

            case 0 : //default case, solve for sensitivity
              sens = sensitivities[least_squares[i].element_number].sensitivity_at_gauss_point[least_squares[i].gauss_point_number];
              break;

            case 1 : // solve for first component of sensitivity
              sens = sensitivities[least_squares[i].element_number].sensitivity_component1_at_gauss_point[least_squares[i].gauss_point_number];
              break;

            case 2 : // solve for second component of sensitivity
      	      sens = sensitivities[least_squares[i].element_number].sensitivity_component2_at_gauss_point[least_squares[i].gauss_point_number];
      	      break;

            }

            // Storing weighted sensitivity.

            B_dash[i] = sens*lsweight;

        }
    }

    // Solve Least Squares without LAPACK
    // step 1 compute A'*A

    int num_pts = number_of_gauss_points ;
    int basis_pts = n;

    std::vector<std::vector<double>> AtA;
    AtA.resize(basis_pts);
    for(int i = 0; i< basis_pts; i++) AtA[i].resize(basis_pts,0.0);

    for(int i = 0; i < basis_pts; i++)
    {
      for(int j = 0; j < basis_pts; j++)
      {

        for(int ii = 0; ii < num_pts; ii++)
        {
          AtA[i][j] += A_dash[ii][i]*A_dash[ii][j];

        }

      }

    }

    // step 2 compute A'*B in vector format
    std::vector<double> AtB(basis_pts,0.0);

    for(int i = 0; i < basis_pts; i++)
    {
      for(int ii = 0; ii < num_pts; ii++)
      {
        AtB[i] += A_dash[ii][i]*B_dash[ii];
      }
    }

    // step 3 compute (A'A)X = A'B using CG


    // initialize every thing here

    int matrix_size = basis_pts;

    std::vector<double> X(matrix_size,0.0);
    std::vector<double> r_inhouse = AtB; // residual

    std::vector<double> p_inhouse = r_inhouse; // direction vector
    double alpha = 0.0;
    double beta = 1.0;

    int max_iter = matrix_size;

    // singularity flag for A'A
    bool is_singular = 0;

    for(int iter = 0; iter < max_iter; iter++)
    {
      std::vector<double> AtAp(matrix_size,0.0);
      AtAp = mat_vec_mult(AtA,p_inhouse); // AtA*p

      double divisor = vec_vec_mult(p_inhouse, AtAp);

      if(std::abs(divisor) == 0){
        is_singular = 1; // this means A'A is singular
        break;
      }

      alpha = vec_vec_mult(r_inhouse,r_inhouse)/divisor ;



      for(int i = 0; i < matrix_size; i++) X[i] += alpha*p_inhouse[i];

      beta = 1.0/vec_vec_mult(r_inhouse,r_inhouse);

      for(int i = 0; i < matrix_size; i++) r_inhouse[i] -= alpha*AtAp[i];

      beta *= vec_vec_mult(r_inhouse,r_inhouse);

      for(int i = 0; i < matrix_size; i++) p_inhouse[i] = r_inhouse[i] + beta*p_inhouse[i];

    }

    if(isnan(X[0]) || is_singular)
    {
      sens = 0.0;
    }
    else
    {
      sens =  X[0];
    }





    return sens;

}

void SensitivityAnalysis :: WriteAverageVonMisesTxt (int datapoint, int num_elem_x, int num_elem_y, std::string file_path, std::string txt_file_name) {

  // Define full file name variables
  ostringstream file_name, num;

  // Create the iteration number
  num.str("");
  num.width(4);
  num.fill('0');
  num << std::right << datapoint;

  file_name.str("");

  if (file_path == "") {

    file_name << txt_file_name << "_" << num.str() << ".txt";

  } else file_name << file_path << "/" << txt_file_name << "_" << num.str() << ".txt";


  FILE *print_file;

  print_file = fopen(file_name.str().c_str(), "w");

  // total number of elements
  int num_elem = num_elem_x * num_elem_y;

  fprintf(print_file, "%i \n", num_elem_x);
  fprintf(print_file, "%i \n", num_elem_y);

  for (int i = 0; i < num_elem; i++) {

    // print each stress value
    fprintf(print_file, "%.16lf \n", strains[i].von_mises_average);

  }

  fclose(print_file);

}

//TODO: write this function
void SensitivityAnalysis :: WriteAverageVonMisesVtk () {
}

void SensitivityAnalysis :: WriteAverageSensitivitiesTxt (int datapoint, int num_elem_x, int num_elem_y, std::string file_path, std::string txt_file_name) {

  // Define full file name variables
  ostringstream file_name, num;

  // Create the iteration number
  num.str("");
  num.width(4);
  num.fill('0');
  num << std::right << datapoint;

  file_name.str("");

  if (file_path == "") {

    file_name << txt_file_name << "_" << num.str() << ".txt";

  } else file_name << file_path << "/" << txt_file_name << "_" << num.str() << ".txt";


  FILE *print_file;

  print_file = fopen(file_name.str().c_str(), "w");

  // total number of elements
  int num_elem = num_elem_x * num_elem_y;

  fprintf(print_file, "%i \n", num_elem_x);
  fprintf(print_file, "%i \n", num_elem_y);

  for (int i = 0; i < num_elem; i++) {

    // print each stress value
    fprintf(print_file, "%.16lf \n", sensitivities[i].sensitivity_average);

  }

  fclose(print_file);

}

// this function multiplies the sparse matrix with a vector
std::vector<double> SensitivityAnalysis ::mat_vec_mult( std::vector<std::vector<double>> &AtA, std::vector<double> &v_in )
{
  // initialze output to a vector of zeros
  std::vector<double> v_out(v_in.size(),0.0);
	//#pragma omp parallel num_threads(4)
  for(int i = 0; i < v_in.size(); i++)
	{
    for(int j = 0; j < v_in.size(); j++)
  	{
      v_out[i] += AtA[i][j]*v_in[j];
    }

	}
  return v_out;
}

// this function multiplies the sparse ve
double SensitivityAnalysis ::vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 )
{
  // initialze output to a vector of zeros
  double result = 0.0;

  for(int i = 0; i < v_in1.size(); i++)  result += v_in1[i]*v_in2[i];

  return result;
}

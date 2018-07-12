
namespace M2DO_LSM {
    
    /*
     ! Structure containing attributes for an individual grid node.
     */
    struct h_Node {
        std::vector<double> sensitivities;  //!< Objective and constraint sensitivities.
    };
    
    /*
     ! Structure containing attributes for an individual grid element.
     */
    struct h_Element {
        std::vector<double> sensitivities;  //!< Objective and constraint sensitivities.
    };
    /*
     Get the index and map of elements and nodes involved in the hole insertion
     */
    int hole_map (M2DO_LSM::Mesh lsmMesh, M2DO_LSM::LevelSet levelSet, double h, double lBand, vector <double>& h_index, vector <bool>& h_elem) {
        
        // Read in mesh variables
        int nNodes = lsmMesh.nNodes;
        int nElementsX = lsmMesh.width;  // Number of elements in X
        int nElementsY = lsmMesh.height; // Number of elements in Y
        int nElements = lsmMesh.nElements;
        vector <double> lsf; // Primary level set function
        
        //double h; // unit length of element ???
        int HC; // parameters for narrow band ???
        double NB = lBand*h;
        
        int nd_count = 0;
        //vector <int> h_index(nNodes); // Identifier of whether this node can have a hole inserted.
        h_index.resize(nNodes); fill(h_index.begin(), h_index.end(), 0);
        for (int inode = 0; inode < nNodes; inode++) {
            
            h_index[inode] = 0;
            
            // Check whether a node is fixed or inside narrow band via its
            //   properties lsmMesh.nodes[theNode].isActive &
            //              lsmMesh.nodes[theNode].isFixed
            //
            // Feasible nodes should be :
            //   a. outside 'narrow band' area;
            //   b. active;
            //   c. not be fixed.

            if ( ( levelSet.signedDistance[inode] >= NB ) && 
                 ( lsmMesh.nodes[inode].isActive || !lsmMesh.nodes[inode].isFixed ) ) {
                nd_count++;
                h_index[inode] = 1;
            }
        }

        // Count the number of elements that can be inserted hole
        int count = 0;
        h_elem.resize(nElements); fill(h_elem.begin(), h_elem.end(), false);
        for (int iel = 0; iel < nElements; iel++) {
            //h_elem[iel] = 0;
            bool isAvailableForHole = false;
            for (int ind = 0; ind < 4; ind++) {
                int inode = lsmMesh.elements[iel].nodes[ind];
                if ( ( levelSet.signedDistance[inode] >= NB ) && 
                     ( lsmMesh.nodes[inode].isActive || !lsmMesh.nodes[inode].isFixed ) ) {
                    isAvailableForHole = true;
                }
            }
            if (isAvailableForHole) { count++; h_elem[iel] = true;}
        }

        cout << "\nNumber of feasible LS nodes are: " << nd_count <<endl;
        cout << "\nNumber of available elements for hole insertion is: " << count << endl;

        return(count);

    }
    
    /*
     Set value of the secondary level set function
     */
    // void get_h_lsf( int nNodes, vector <double>& h_index, vector <M2DO_LSM::h_Node> h_Nsens,
    //                vector <double> gammas, vector <double>& h_lsf ) {
        
    //     // h_lsf.resize(nNodes); fill( h_lsf.begin(), h_lsf.end(), 0.0 );
        
    //     for ( int inode = 0; inode < nNodes; inode++ ) {
    //         // h_lsf[inode] = 1.0;
            
    //         if ( h_index[inode] ) {
                
    //             for (int j = 0; j < h_Nsens[inode].sensitivities.size(); j++ ) {
    //                 h_lsf[inode] -= gammas[j] * ( 1 - h_Nsens[inode].sensitivities[j] );
    //             }
    //             if (h_lsf[inode]<0) { 
    //                 cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    //                 cout << "\n  The h_lsf of node " << inode << "\t is " << h_lsf[inode] << endl;
    //                 cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl; 
    //             }
    //         }
            
    //     }
        
    // }

    void get_h_lsf( int nNodes, vector <double>& h_index, vector <M2DO_LSM::h_Node> h_Nsens,
                   vector <double> gammas, vector <double>& h_lsf ) {
        
        // h_lsf.resize(nNodes); fill( h_lsf.begin(), h_lsf.end(), 0.0 );
        
        for ( int inode = 0; inode < nNodes; inode++ ) {
            // h_lsf[inode] = 1.0;
            
            if ( h_index[inode] ) {
                
                for (int j = 0; j < h_Nsens[inode].sensitivities.size(); j++ ) {
                    h_lsf[inode] -= gammas[j] * h_Nsens[inode].sensitivities[j];
                }
                // if (h_lsf[inode]<0) { 
                //     cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                //     cout << "\n  The h_lsf of node " << inode << "\t is " << h_lsf[inode] << endl;
                //     cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl; 
                // }
            }
            
        }
        
    }

    void get_h_lsf2( int nNodes, vector <double>& h_index, vector <M2DO_LSM::h_Node> h_Nsens,
                   vector <double> gammas, vector <double>& h_lsf ) {
        
        // h_lsf.resize(nNodes); fill( h_lsf.begin(), h_lsf.end(), 0.0 );
        
        for ( int inode = 0; inode < nNodes; inode++ ) {
            h_lsf[inode] = 1.0;
            
            if ( h_index[inode] ) {
                
                for (int j = 0; j < h_Nsens[inode].sensitivities.size(); j++ ) {
                    h_lsf[inode] -= gammas[j] * (1 - h_Nsens[inode].sensitivities[j]);
                }
                // if (h_lsf[inode]<0) { 
                //     cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
                //     cout << "\n  The h_lsf of node " << inode << "\t is " << h_lsf[inode] << endl;
                //     cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl; 
                // }
            }
            
        }
        
    }
    
    /*
    Get area of elements availabe to insert hole
    */
    double get_h_area( M2DO_LSM::Mesh lsmMesh, vector<bool>& h_elem ) {
    
        double Atotal = 0.0;

        for ( int iel = 0; iel < lsmMesh.nElements; iel++ ) { 

            if ( h_elem[iel] ) { // available element for hole insertion
                // Add to the total area
                Atotal += lsmMesh.elements[iel].area;
            }
        }
        cout << "\nAvailable area for hole insertion is:\t " << Atotal << endl;
        return(Atotal);

    }

    /*
        Compute area related to holes
    */

    // double computeAreaFractions_hole(){
    //     //
    // }    

    /*
    Function to find minimal value of h_lsf, thus revealing if signed distance 
    function has been changed beyond CFL condition
    */
    double find_min_h_lsf(vector <double>& h_lsf) {

        double Atemp = 1000;
        int nNodes = h_lsf.size();

        for ( int inode = 0; inode < nNodes; inode++ ) {
            Atemp = h_lsf[inode] < Atemp ? h_lsf[inode] : Atemp;
        }

        return Atemp;

    }

    /*
    Function to do Finite Difference analysis on minimum value of h_lsf
    */
    double hole_FD_h_lsf(vector <double>& h_index, 
                         vector <M2DO_LSM::h_Node> h_Nsens,
                         vector <double>& h_lsf,
                         vector <double> gammas,
                         double DA,
                         double TargetV,
                         double perb) {

        double Atemp, DAperb, FDiff;
        int nGammas = gammas.size();
        int nNodes = h_lsf.size();

        // Get perturbed gamma
        for ( int i = 0; i < nGammas; i++ ) {
            gammas[i] = gammas[i] + perb;
        }

        // Get hole lsf and total area of perturbed gamma
        get_h_lsf2( nNodes, h_index, h_Nsens, gammas, h_lsf );
        Atemp = find_min_h_lsf( h_lsf );

        // Get change in the objective
        DAperb = Atemp - TargetV;

        FDiff = ( DA - DAperb ) / perb;

        return(FDiff);
    }

    /*
     Function to initialise the hole level set function and the maximum and
     minimum values of gamma (side limits)
     */
    
    void initialise_hole_lsf(M2DO_LSM::Mesh lsmMesh, int h_count, double holeCFL,
                             M2DO_LSM::LevelSet levelSet,
                             double moveLimit, 
                             vector <double>& h_index,
                             vector <bool>& h_elem,
                             vector <M2DO_LSM::h_Node> h_Nsens,
                             vector <double>& h_lsf,
                             vector <double> lambdas) {
        
        // read mesh data
        double h; // unit length of element in level set mesh (?)
        int nNodes = lsmMesh.nNodes; // Number of nodes in level set mesh
        int nElements = lsmMesh.nElements; // Number of elements in level set mesh
        // vector <double> h_lsf(nNodes); // Hole level set function
        
        // vector <double> h_Nsens; // Node sensitivity
        // vector <double> h_Esens; // Element sensitivity
        
        int nGammas = lambdas.size(); // Number of constraints + 1
        vector <double> gammas(nGammas); // Weight factors for the linear combination of sensitivities
        vector <double> h_area;  //  Store change in element area caused by holes
        vector <double> h_gMin(nGammas);  // Minimum hole weighting
        vector <double> h_gMax(nGammas);  // Maximum hole weighting
        
        /*
         int *h_index = malloc(inMesh.NumNodes*sizeof(int)); //Identifier of weather this node can have a hole inserted
         int *h_EmapX = malloc(inMesh.NumElem*sizeof(int));  //Element mapping in x direction
         int *h_EmapY = malloc(inMesh.NumElem*sizeof(int));  //Element mapping in y direction
         int *h_posN = malloc((1+lsprob.num)*sizeof(int));   //Node sensitivity mapping
         int *h_posE = malloc((1+lsprob.num)*sizeof(int));   //Element sensitivity mapping
         double *h_Nsens = malloc((1+lsprob.num)*inMesh.NumNodes*sizeof(double)); //Node sensitivity
         double *h_Esens = malloc((1+lsprob.num)*inMesh.NumElem*sizeof(double)); //Element sensitivity
         double *h_lsf = malloc(inMesh.NumNodes*sizeof(double)); //hole level set funtion
         double *h_area = malloc(inMesh.NumElem*sizeof(double)); //store change in element area caused by holes
         double *h_gMin = malloc(1*sizeof(double));  //Minimumum Hole weighting
         double *h_gMax = malloc(1*sizeof(double));  //Maximumum Hole weighting
         int *Reint = malloc(1*sizeof(int));         //Flag to let us know if we need to reinitailise the model
         */
        
        double minlsf = 0.0, max_dh_lsh;
        double dh_lsf, temp, Atemp;
        
        /*
         Find minimum of h_lsf or maximum of change of h_lsf,
           and in turm get minmum gamma values ??? [for setting side limits for gamma]
         */
        for (int i = 0; i < nGammas; i++ ) { gammas[i] = 1.0; } 

        int h_nodes_count = 0;
        for (int inode = 0; inode < nNodes; inode++) {
            
            h_lsf[inode] = 1.0;
            
            if (h_index[inode]) {

                h_nodes_count++; 
                
                dh_lsf = 0;
                
                for (int j = 0; j < nGammas; j++) {

                    dh_lsf += gammas[j] * h_Nsens[inode].sensitivities[j];
                }
                if (h_nodes_count == 1) {
                    max_dh_lsh = dh_lsf; 
                } else {
                    max_dh_lsh = ( max_dh_lsh > dh_lsf ) ? max_dh_lsh : dh_lsf ;
                }
                
                //cout<<"Node " << inode << " minlsf is " << minlsf << " delta lsf is " << temp << " h_lsf " << h_lsf[inode] << endl;
            }
        }

        h_gMin[0] = (1.0/max_dh_lsh)*1.00;
        cout<<"\nMinimal weight factor should be: " << h_gMin[0] << endl;
        
        for (int i = 0; i < nGammas; i++ ) { gammas[i] = h_gMin[0]; } 
        get_h_lsf2( nNodes, h_index, h_Nsens, gammas, h_lsf );
        
        // compute the sum of availabe area for hole insertion
        // Atemp = get_h_area( lsmMesh, h_elem );
        
        // Now use Newton's method to find maximum of gammer or lambda
        double Gtemp = 1.0001 * h_gMin[0];
        double perb = 0.0001 * h_gMin[0];
        double Gabove = 10000 * Gtemp;
        double Gbelow = 0.0;
        double gamma_old = Gtemp;
        
        double Atotal = h * h * h_count;
        double TargetV = -1.0 * holeCFL * h; // TargetV = Atotal * ( 1 - holeCFL );
        double DA, DAP;
        double FDiff;
        
        int loop = 0;
        do {
            gamma_old = Gtemp;
            get_h_lsf2( nNodes, h_index, h_Nsens, gammas, h_lsf );
            Atemp = find_min_h_lsf( h_lsf );
        
            // See if removed area Atotal is the matches the targeted maximum (within range)
            DA = Atemp - TargetV;
            DAP = Atemp / Atotal;
        
            // Update overshoot variables
            if ( loop > 2 ) {
                if ( DA > 0 ) { 
                    Gbelow = ( Gtemp > Gbelow ) ? Gtemp : Gbelow; 
                }
                else if ( DA < 0 ) {
                    Gabove = ( Gtemp < Gabove ) ? Gtemp : Gabove;
                }
            }
            
            // Make allowances for initial gradient being so low
            if ( Atotal == 0 ) { Gabove = Gtemp - perb; }
        
            // some words here
            if ( fabs( DA ) < 0.01 ) {
                loop = 501;
            }
            else {
                // Get gradient by Finite Difference
                FDiff = hole_FD_h_lsf(h_index, h_Nsens, h_lsf, gammas, DA, TargetV, perb);
        
                // Cover for no gradient due to too small Gtemp (gammer temp)
                if (FDiff == 0) {
                    Gtemp += perb;
                }
                else {
                    Gtemp += ( DA / FDiff );
                }
        
                // If outsie range then biset the range for the new value.
                if ( ( Gtemp <= Gbelow) || ( Gtemp >= Gabove ) || ( Atemp = 0.0 ) ) {
                    Gtemp = ( Gabove + Gbelow ) / 2.0;
                }
        
                // Update Gtemp (gammer temp)
                for ( int j = 0; j < nGammas; j++ ) {
                    gammas[j] = Gtemp;
                }
        
                loop++;
            }
        } while( loop < 500 );
        
        if ( loop == 500 ) {
            printf("\nWARNING HOLE CFL DID NOT CONVERGE, HOLES MAY BE TOO LARRGE");
        } 

        // Assign value for upper bound of gamma
        h_gMax[0] = Gtemp;
 
        printf( "\nHole gamma side limits: %f and %f \n.", h_gMin[0], h_gMax[0] );
        get_h_lsf2( nNodes, h_index, h_Nsens, gammas, h_lsf );
        // // For plotting set up the output to be gammer max
        // for ( int j = 0; j < num_gam; j++ ) { Gam[j] = h_gMax[0]; }
        // get_h_lsf( nNodes, h_index, h_Nsens, lambdas, h_lsf );
        // Atemp = M2DO_LSM::Boundary::computeAreaFractions();
    }
    
    
    // /*
    //     Conduct the Finite Difference method to calculate gradient
    // */
    // double hole_FD_SG( int nNodes, vector <double>& h_index,
    //                   vector <M2DO_LSM::h_Node> h_Nsens, vector <double> lambdas,
    //                   vector <double>& h_lsf, double perb, double DA, doulbe TargetV ) {
    
    //     double Atemp, DAperb, FDiff;
    //     int num_gam = lambdas.size();
    
    //     // Get perturbed gammer
    //     for ( int i = 0; i < num_gam; i++ ) {
    //         Gam[i] = Gam[i] + perb;
    //     }
    
    //     // Get hole lsf and total area of perturbed gammer
    //     get_h_lsf( nNodes, h_index, h_Nsens, lambdas, h_lsf );
    //     Atemp = M2DO_LSM::Boundary::computeAreaFractions();
    
    //     // Get change in the objective
    //     DAperb = Atemp - TargetV;
    //     FDiff = ( DA - DAperb ) / perb;
    
    //     return(FDiff);
    
    // }
    
    
    // // Calculate sensitivities using least squares of integratiin points for AFG
    // // method
    
    // void sens_hole();
    
    // // Function to set velocity (or move dist) using SLP (???)
    // // - filter method with jole insertion
    
    // void SLPsubSol4_hole();
    /*
    void reinitialise_hole(){

        int nNodes = lsmMesh.nNodes;

        // 
        for (int inode = 0; inode < nNodes; inode++) {
            lsf_temp[inode] = 0.0;
        }

        // Ensure all nodes on domain boundary have lsf <=0.0

        // 
        // Initalise Known set
        for(int inode = 0; inode < nNodes; inode++) {
            // If node is fixed or on boundary, then retain the lsf
            if(fabs(lsf[inode]) < tol) {
                known[i] = true;
                lsf_temp[inode] = lsf[inode]; // or 0.0 ??
            }
        }

        // 
        for (int sign = 1; sign > -2; sign-=2) {
            ftemp = h * 1000.0;
            // Initialise trial set
            for(inode = 0; inode < nNodes; inode++) {
                    
                // If current node not in known set and is the correct side of the boundary
                if(!known[inode] && (sign * lsf[inode]) > 0.0) {
                    NS = ftemp * lsf[inode];
                    // If any neighbouring node is on opposite side of the boundary, or on boundary, then it is a trial node
                    for (int i = 0; i < 4; i++) {
                        int neighbour_node;
                        neighbour_node = lsmMesh.nNodes[inode].neighbours[i];
                        if (lsf[neighbour_node] * NS < 0.0) { trial[inode] = true; }
                    }
                }
            }

            // Calculate lsf_temp for all current trial nodes
            // ???LocalInt

            // Update trial set by considering neighbours of nodes in known set
            for (int inode = 0; inode < nNodes; inode++ ) {
                if ( !known[inode] && (sign * lsf[inode]>0.0) ) {
                    // If any neighbouring node not in known set and is the current side of the boundary
                    for (int i = 0; i < 4; i++) {
                            int neighbour_node;
                            neighbour_node = lsmMesh.nNodes[inode].neighbours[i];
                            if (known[neighbour_node]) { trial[inode] = true; }
                    }
                }
            }

            // Do until trial set is empty
            int flag, ncount;
            vector <double> lsf_trial(nNodes);
            vector <int> xind;
            do {
                flag = 0;
                // for (int inode = 0; inode < nNodes; inode++) {
                //     flag = 0;
                //     if (trial[inode]) {
                //         flag = 1;
                //         if ( fabs(lsf_trial[inode]) < 0.0 ) {
                //             Af = 0.0; Bf = 0.0; Cf = 0.0;
                //             ftemp = 0.0; // initial


                //         }

                //     }
                // }

                if (flag == 0) { break; } // if there are no trial nodes, then end the loop

                ncount = 0; // re-initialise
                // check to see which nodes need updating this iteration
                for (int inode = 0; inode < nNodes; inode++) {

                    if ( lsf_trial[inode] < 0.0 ) {
                        printf("\nERROR IN ReInt:");
                    }

                    // If a trial node has a trial lsf value <= minimum then update
                    if ( trial[inode] && (fabs(lsf_trial[inode] - lsf_min) < tol ) ) {
                        lsf_temp[inode] = sign * lsf_trial[inode];
                        // Move node to known set
                        xind.pushback(inode);
                    }
                }

                // Update knonw set for next iteration
                // .1
                for (int i = 0; i<ncount; i++) {
                    trial[xind[i]] = false;
                    knonw[xind[i]] = true;
                }

                // .2
                for (int i = 0; i<ncount; i++) {
                    int inode_temp = xind[i];
                    for (int j = 0; j < 4; j++){
                        int neighbour_node = lsmMesh.nodes[inode_temp].neighbours[j];
                        if (!known[neighbour_node] && (sign * lsf[neighbour_node] > 0.0) ) {
                            trial[neighbour_node] = true;
                            lsf_trial[neighbour_node] = 0.0;
                        }
                    }
                }

                if (ncount == 0) {
                    printf("\nERROR! ncount=0 in ReInt! Aborting");
                }

            } while ( ncount > 0 )
        }

        
    }*/

    
}


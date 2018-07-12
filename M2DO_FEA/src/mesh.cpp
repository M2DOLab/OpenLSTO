#include "quadrature.h"
#include "linear_shape_function.h"
#include "node.h"
#include "solid_element.h"
#include "solid_material.h"
#include "mesh.h"

using namespace M2DO_FEA ;

Mesh :: Mesh () {

	is_structured = false ; // By default.

}

Mesh :: Mesh (int spacedim) : spacedim (spacedim) {

	is_structured = false ; // By default.

}

void Mesh :: Print () {

	cout << "Mesh (" ;

	for (int i = 0 ; i < nodes.size() ; ++i) {

		if (i > 0) {

			cout << ", " ;

		}

		nodes[i].Print() ;
	}

	cout << ")" ;

}

void Mesh :: MeshSolidHyperRectangle (vector<int> nel, MatrixXd mesh_box, int element_order, bool time_it) {

	auto t_start = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "\nMeshing solid hyper-rectangle ... " << flush ;

	}

	/*
		First mesh the natural hyper-rectangle. We can utilise a sequence
		similar to that in quadrature.UpdateEtaCounter() for this.
	*/

	int num_nodes = 1, num_elements = 1 ;

	for (int i = 0 ; i < nel.size() ; ++i) {

		num_nodes    *= nel[i]+1 ;
		num_elements *= nel[i] ;

	}

  	nodes.reserve (num_nodes) ;
  	solid_elements.reserve (num_elements) ;

	/*
		Each of these linear elements has pow (2, spacedim) nodes;
		each node has dim degrees of freedom. Using this we can
		calculate num_entries (needed for study step).
	*/

	// num_entries += num_elements * pow((pow(2, spacedim) * dim), 2) ;
	// n_dof 	     = num_nodes * dim ;

	vector<int> eta_count (spacedim, 0) ;

	for (int i = 0 ; i < num_nodes ; ++i) {

		Node node (spacedim) ;
		node.id = i ;

		for (int l = 0 ; l < spacedim ; ++l) {

			node.coordinates[l] = -1 + eta_count[l] * (2.0 / nel[l]) ;

		}

		nodes.push_back(node) ;

		// Update eta_count:

		eta_count[0] += 1 ;

		if (eta_count[0] > nel[0]) {

			eta_count[0] = 0 ;

			for (int l = 1 ; l < spacedim ; ++l) {

				eta_count[l] += 1 ;

				if (eta_count[l] <= nel[l]) {
					break ;
				}

				else {
					eta_count[l] = 0 ;
				}

			}

		}

	}

	/*
		We create the elements on the natural hyper_rectangle,
		as it is easy to visualize. The node locations will change
		below, but that doesn't affect the elements; they need
		only know which nodes belong to them.
	*/

	LinearShapeFunction linear_shape_function (spacedim, 1) ;

	vector<int> node_ids (pow(2, spacedim), 0) ;
	vector<int> nel_count (spacedim, 0) ;
	vector<int> nel_mult (spacedim, 1) ;
	vector<double> eta (spacedim, 0) ;
	int eta_int, start_id ;

	for (int l = 1 ; l < spacedim ; ++l) {

		nel_mult[l] *= nel_mult[l-1] * (nel[l-1]+1) ;

	}

	SolidElement element (spacedim, element_order, *this) ;

	for (int i = 0 ; i < num_elements ; ++i) {

		start_id = 0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			start_id += nel_count[l] * nel_mult[l] ;

		}

	// This just fills the node_ids vector with start_id:
		fill(node_ids.begin(), node_ids.end(), start_id) ;

		for (int j = 0 ; j < pow(2, spacedim) ; ++j) {

			eta = linear_shape_function.GetEta(j) ;

			// Change the -1 values to zeros. Also, eta
			// comes as doubles so change to int, then
			// multiply by nel_mult:

			for (int k = 0 ; k < spacedim ; ++k) {

				eta_int = (eta[k] < 0) ? 0 : 1 ;
				node_ids[j] += eta_int * nel_mult[k] ;

			}

			element.node_ids[j] = node_ids[j] ;

		}

		// Add the element to the mesh:
		solid_elements.push_back(element) ;

		// Update nel_count:

		nel_count[0] += 1 ;

		if (nel_count[0] > nel[0]-1) {

			nel_count[0] = 0 ;

			for (int l = 1 ; l < spacedim ; ++l) {

				nel_count[l] += 1 ;

				if (nel_count[l] < nel[l]) {
					break ;
				}

				else {
					nel_count[l] = 0 ;
				}

			}

		}


	}


	/*
		Now we deform the natural mesh to conform to the physical
		geometry using linear shape functions:

		x = sum(N_i * x_i) etc. where in this case x_i are the
		coordinates of the corners of the hyper_rectangle to mesh.
	*/

	VectorXd shape_value_vec ;

	for (int i = 0 ; i < nodes.size() ; ++i) {

		// Shape function values at natural coordinate:
		shape_value_vec = linear_shape_function.GetShapeFunctionValuesVector(nodes[i].coordinates) ;

		for (int j = 0 ; j < spacedim ; ++j) {

			nodes[i].coordinates[j] = mesh_box.col(j).dot(shape_value_vec)  ;

		}

	}

	/*
		Now compute centroids:
	*/

	ComputeCentroids () ;

	/*
		Compute total time:
	*/

	auto t_end = chrono::high_resolution_clock::now() ;

	if (time_it) {

		cout << "Done. Time elapsed = " << chrono::duration<double>(t_end-t_start).count() << "\n" << flush ;

	}

}

void Mesh :: AssignDof () {

	n_dof = 0 ;

	/*
		Solid elements: these have spacedim displacement dofs per node,
		and there are pow(2, spacedim) nodes per element.
	*/

	for (auto && element : solid_elements) {

		element.dof = vector<int> (spacedim * pow (2, spacedim), -1) ;

		for (int i = 0 ; i < element.node_ids.size() ; ++i) {

			auto && node = nodes[element.node_ids[i]] ;

			/*
				Check if this node has already assigned a dof number;
				If not, assign one:
			*/

			for (int j = 0 ; j < spacedim ; ++j) {

				if (node.dof[j] >= 0) {

					element.dof [i*spacedim + j] = node.dof [j] ;

				}

				else {

					element.dof [i*spacedim + j] = n_dof ;
					node.dof [j] = n_dof ;
					n_dof += 1 ;

				}

			} // for j (spacedim).

		} // for i (element.node_ids.size()).

	} // for element : solid_elements.

}

int Mesh :: n_entries () {

	int count = 0 ;

	count += solid_elements.size() * pow(spacedim * pow(2, spacedim), 2) ;

	// etc.

	return count ;

}

// int Mesh :: nc_entries () {

// 	int count = 0 ;

// 	count += poisson_elements.size() * pow(pow(2, spacedim), 2) ;

// 	// etc.

// 	return count ;

// }

vector<int> Mesh :: GetNodesByCoordinates (vector<double> coord, vector<double> tol) {

	vector<int> selected_nodes ;
	bool inside ;

	for (int i = 0 ; i < nodes.size() ; ++i) {

		inside = true ;

		for (int j = 0 ; j < spacedim ; ++j) {

			if ( (abs(nodes[i].coordinates[j] - coord[j]) > tol[j]) ) {

				inside = false ;

			}

		}

		if ( inside ) {
			selected_nodes.push_back (nodes[i].id) ;
		}

	}

	return selected_nodes ;

}

vector<int> Mesh :: dof (int node_id) {

	int dof_count = 0 ;
	vector<int> dof_vec (6, -1) ;

	for (auto && d : nodes[node_id].dof) {

		if (d >= 0) {

			dof_vec[dof_count] = d ;
			dof_count += 1 ;

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

vector<int> Mesh :: dof (vector<int> node_ids) {

	int dof_count = 0 ;
	vector<int> dof_vec (6 * node_ids.size(), -1) ;

	for (auto && id : node_ids) {

		for (auto && d : nodes[id].dof) {

			if (d >= 0) {

				dof_vec[dof_count] = d ;
				dof_count += 1 ;

			}

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

vector<int> Mesh :: dof (vector<int> node_ids, vector<int> components) {

	int dof_count = 0 ;
	vector<int> dof_vec (components.size() * node_ids.size(), -1) ;

	for (auto && id : node_ids) {

		for (auto && c : components) {

			int d = nodes[id].dof[c] ;

			if (d >= 0) {

				dof_vec[dof_count] = d ;
				dof_count += 1 ;

			}

		}

	}

	dof_vec.resize (dof_count) ;

	return dof_vec ;

}

void Mesh :: ComputeCentroids () {

	vector<double> eta (spacedim, 0.0);

	for (int i = 0 ; i < solid_elements.size() ; i++) {

		auto && element = solid_elements[i] ;

		element.centroid.resize (spacedim) ;

		VectorXd xc = element.NaturalToPhysicalCoordinates (eta) ;

		for (int k = 0 ; k < spacedim ; k++) {

			element.centroid[k] = xc[k] ;

		}

	}

}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Mesh :: nodeConnectivity() {

node_connectivity.resize(nodes.size());

for (int i = 0 ; i < nodes.size(); i++){


		for (int j = 0 ; j < solid_elements.size(); j++){

			auto && element = solid_elements[j] ;


          for (int k = 0; k < element.node_ids.size() ; k++){


			if (element.node_ids[k] == i) node_connectivity[i].push_back(j);


			}

		  }

		}
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Mesh :: elementNeighbours() {

element_neighbours.resize(solid_elements.size());


for (int i = 0 ; i < solid_elements.size(); i++){

	auto && element = solid_elements[i] ;


	// nodes_of_element_i.resize(element.node_ids.size())	;

    vector<int> nodes_of_element_i = element.node_ids;
	//cout<<endl<<"element "<<i<<endl;

	for (int j = 0; j < solid_elements.size(); j++){



		auto && element1 = solid_elements[j] ;

		// nodes_of_element_j.resize(element1.node_ids.size());

        vector<int> nodes_of_element_j = element1.node_ids;
//cout<<"HI"<<endl;

 		int r = 0;

        if (i != j){



			for (int k = 0 ; k < element.node_ids.size(); k++){

				for (int h = 0; h < element1.node_ids.size(); h++){



				//cout<<"Hello"<<endl;
				if (nodes_of_element_i[k] == nodes_of_element_j[h]) {

					//cout<<" shares with neighbour "<<j<<" the nodes: "<<endl;

					//cout<<nodes_of_element_j[h]<<endl; //this will find all the nodes that the element shares with its surrounding elements
					//including the diagonal elements (i.e not only side neighbours but also corner neighbours)

					r += 1;
				//cout<<r<<endl;
				}


			}


			}



				}
			if (r == 2) {
					element_neighbours[i].push_back(j);


				//cout<<endl<<" neighbours: "<<j<<endl;

			}
			//cout<<r<<endl;

			}

			//


	}

	//for (int i = 1; element_neighbours.size(); i++){

	//	cout<<"Hello"<<endl;
	//}



}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


// Write node coordinates and element connectivity.
void Mesh :: WriteMeshTXT() {

	// NODES

    // Defining file name variables.
    ostringstream fileName;

    // Creating file name.
    fileName.str("");
    fileName << "Output/nodes.txt";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    for (int i = 0; i < nodes.size(); i++)
    {
    	for (int j = 0; j < spacedim; j++)
    	{
			fprintf(pFile, "%lf \t", nodes[i].coordinates[j]);
    	}
    	fprintf(pFile, "\n");
    }
    fclose(pFile);

    // ELEMENTS

    // Defining file name variables.
    ostringstream fileName_el;

    // Creating file name.
    fileName_el.str("");
    fileName_el << "Output/elements.txt";

    FILE *pFile_el;

    pFile_el = fopen(fileName_el.str().c_str(), "w");

    for (int i = 0; i < solid_elements.size(); i++)
    {
    	fprintf(pFile, "%li \t", element_type[i]);
    	for (int j = 0; j < solid_elements[i].node_ids.size(); j++)
    	{
			fprintf(pFile, "%li \t", solid_elements[i].node_ids[j]);
    	}
    	fprintf(pFile, "\n");
    }
    fclose(pFile);
}

// Write area fractions.
void Mesh :: WriteAreaFractionsTXT(int datapoint, int nelx, int nely) {

    // Defining file name variables.
    ostringstream fileName, num;

    // Creating file name.
    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    fileName << "Output/area_fractions_" << num.str() << ".txt";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // total number of elements.
    int nel = nelx*nely;

    fprintf(pFile, "%i \n", nelx);
    fprintf(pFile, "%i \n", nely);

    for (int i = 0; i < nel; i++)
    {
        // Printing each stress value.
        fprintf(pFile, "%lf \n", solid_elements[i].area_fraction);

    }

    fclose(pFile);
}

// Write area fractions.
void Mesh :: WriteElementTypesTXT(int datapoint, int nelx, int nely) {

    // Defining file name variables.
    ostringstream fileName, num;

    // Creating file name.
    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    fileName << "Output/element_types_" << num.str() << ".txt";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // total number of elements.
    int nel = nelx*nely;

    fprintf(pFile, "%i \n", nelx);
    fprintf(pFile, "%i \n", nely);

    for (int i = 0; i < element_type.size(); i++)
    {
        // Printing each stress value.
        fprintf(pFile, "%i \n", element_type[i]);
    }

    fclose(pFile);
}

void Mesh :: saveNodalPropertiesVTK(const unsigned int& datapoint, int nelx, int nely) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    fileName << "nodal-properties_" << num.str() << ".vtk";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // Set up ParaView header information.
    fprintf(pFile, "# vtk DataFile Version 3.0\n");
    fprintf(pFile, "Para0\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + nelx, 1 + nely, 1);
    fprintf(pFile, "X_COORDINATES %d int\n", 1 + nelx);
    for (unsigned int i=0;i<=nelx;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + nely);
    for (unsigned int i=0;i<=nely;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");
    fprintf(pFile, "POINT_DATA %d\n", nodes.size());

    // Write the nodal signed distance to file.
    fprintf(pFile, "SCALARS distance float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<nodes.size();i++)
        fprintf(pFile, "%lf\n", nodes[i].property);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}

void Mesh :: saveNodalValuesVTK(const unsigned int& datapoint, int nelx, int nely, vector<double> nodal_values) const
{
    std::ostringstream fileName, num;

    num.str("");
    num.width(4);
    num.fill('0');
    num << std::right << datapoint;

    fileName.str("");
    fileName << "nodal-values_" << num.str() << ".vtk";

    FILE *pFile;

    pFile = fopen(fileName.str().c_str(), "w");

    // Set up ParaView header information.
    fprintf(pFile, "# vtk DataFile Version 3.0\n");
    fprintf(pFile, "Para0\n");
    fprintf(pFile, "ASCII\n");
    fprintf(pFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(pFile, "DIMENSIONS %d %d %d\n", 1 + nelx, 1 + nely, 1);
    fprintf(pFile, "X_COORDINATES %d int\n", 1 + nelx);
    for (unsigned int i=0;i<=nelx;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nY_COORDINATES %d int\n", 1 + nely);
    for (unsigned int i=0;i<=nely;i++) {fprintf(pFile, "%d ", i);}
    fprintf(pFile, "\nZ_COORDINATES 1 int\n0\n\n");
    fprintf(pFile, "POINT_DATA %d\n", nodes.size());

    // Write the nodal signed distance to file.
    fprintf(pFile, "SCALARS distance float 1\n");
    fprintf(pFile, "LOOKUP_TABLE default\n");
    for (unsigned int i=0;i<nodes.size();i++)
        fprintf(pFile, "%lf\n", nodal_values[i]);

    fclose(pFile);

    return;

error:
    exit(EXIT_FAILURE);
}



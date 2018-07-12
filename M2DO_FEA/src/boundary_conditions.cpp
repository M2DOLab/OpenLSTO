#include "boundary_conditions.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>

using namespace M2DO_FEA ;

// Function to zip vectors, i.e., store their elements by pairs.
template <typename A, typename B>
void zip(
    const std::vector<A> &a, 
    const std::vector<B> &b, 
    std::vector<std::pair<A,B>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

// Function to unzip vectors, i.e., rewriting the original values.
template <typename A, typename B>
void unzip(
    const std::vector<std::pair<A, B>> &zipped, 
    std::vector<A> &a, 
    std::vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

DirichletBoundaryConditions :: DirichletBoundaryConditions () {

	//

}

DirichletBoundaryConditions :: DirichletBoundaryConditions (vector<int> & dof_in, vector<double> & amplitude_in, int mesh_n_dof_in) : mesh_n_dof (mesh_n_dof_in) {

	if (dof_in.size() != amplitude_in.size())
	{
		cout << "Vectors for DOF's and amplitudes are not the same size!\n";
	}

	// Zip dof_in and amplitude_in vectors together.
    vector<pair<int,double>> zipped;
    zip(dof_in, amplitude_in, zipped);

    // Sort the vector of pairs.
    sort(begin(zipped), end(zipped), 
        [&](const auto& a, const auto& b)
        {
            return a.first < b.first;
        });

    // Write the sorted pairs back to the original vectors.
    unzip(zipped, dof_in, amplitude_in);

    // Store boundary conditions.
    dof = dof_in ;
    amplitude = amplitude_in ;
	MapReducedDofToDof () ; // Compute reduced (active) DOF's.


}

void DirichletBoundaryConditions :: Print () {

	cout << "HomogeneousDirichletBoundaryConditions ( " ;
	cout << "dof = " ;
	// print_vector(dof) ;
	cout << " )" ;

}

void DirichletBoundaryConditions :: MapReducedDofToDof () {

	vector<int> all_dof (mesh_n_dof) ;

	for (int i = 0 ; i < mesh_n_dof ; ++i) {

		all_dof[i] = i ;

	}

	// Identifying dof's with homogeneous boundary conditions (= 0).
	for (int i = 0 ; i < dof.size() ; ++i)
	{
		if (amplitude[i] == 0)
		{
			dof_zeros.push_back(dof[i]);
		}
	}

	set_difference(all_dof.begin(), all_dof.end(), dof_zeros.begin(), dof_zeros.end(), inserter (reduced_dof_to_dof_map, reduced_dof_to_dof_map.begin())) ;

	dof_to_reduced_dof_map.resize (mesh_n_dof) ;

	for (int i = 0 ; i < mesh_n_dof ; ++i) {

		dof_to_reduced_dof_map[i] = -1 ;

	}
	
	for (int i = 0 ; i < reduced_dof_to_dof_map.size() ; ++i) {

		dof_to_reduced_dof_map [reduced_dof_to_dof_map[i]] = i ;

	}

}

PointValues :: PointValues (vector<int> & dof_in, vector<double> & values_in) {

	dof = dof_in ;
	values = values_in ;

	// sort(selected_dof.begin(), selected_dof.end()) ;

}

void PointValues :: Print () {

	cout << "PointValues ( " ;
	cout << "dof = " ;
	// print_vector(dof) ;
	cout << ", values = " ;
	// print_vector(values) ;
	cout << " )" ;

}




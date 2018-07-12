#include "node.h" 

using namespace M2DO_FEA ;

Node :: Node (int spacedim) : spacedim (spacedim), coordinates (spacedim, 0) {
	
	// Constructor puts the node at origin by default.

	// I set the id = -1 to cause failure if this is not set
	// properly at the meshing step.
	id = -1 ;
	dof = vector<int> (7, -1) ;

}

Node :: Node (int spacedim, int id, vector<double> coordinates) : spacedim (spacedim), id (id), coordinates (coordinates) {
	
	dof = vector<int> (7, -1) ;

}

void Node :: Print () {

	cout << "Node (" ;

	for (int i = 0 ; i < coordinates.size() ; ++i) {
		if (i > 0) {
			std::cout << ", " ;
		}

		cout << coordinates[i] ;
	}

	cout << ")" ;

}

vector<double> Node :: ReturnFirstNCoordinates (int n) {

	vector<double> first_coord (n, 0.0) ;

	for (int i = 0 ; i < n ; ++i) {
		first_coord[i] = coordinates[i] ;
	}

	return first_coord ;

}

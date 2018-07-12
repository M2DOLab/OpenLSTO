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

#include <../vendor/eigen3/Eigen/Dense>
#include <../vendor/eigen3/Eigen/Sparse>

template <class T>
void print_vector (vector<T> v) {
	
	cout << "[" ;
	
	for (int i = 0 ; i < v.size() ; ++i) {
		
		if (i > 0) {
			
			cout << "    " ;

		}

		cout << v[i] ;
	}

	cout << "]" ;

}

VectorXd vector_cross (VectorXd & u, VectorXd & v) {

	assert (u.size() == 3) ;
	assert (v.size() == 3) ;

	VectorXd w = VectorXd::Zero(u.size()) ;

	for (int i = 0 ; i < w.size() ; ++i) {

		for (int j = 0 ; j < w.size() ; ++j) {

			if (i == j) {
				// Do nothing.
			}

			else if (i == 0 && j == 1) {
				w(2) += u(i) * v(j) ;
			}

			else if (i == 0 && j == 2) {
				w(1) -= u(i) * v(j) ;
			}

			else if (i == 1 && j == 0) {
				w(2) -= u(i) * v(j) ;
			}

			else if (i == 1 && j == 2) {
				w(0) += u(i) * v(j) ;
			}

			else if (i == 2 && j == 0) {
				w(1) += u(i) * v(j) ;
			}

			else if (i == 2 && j == 1) {
				w(0) -= u(i) * v(j) ;
			}

		}

	}

	return w ;
}

Vector3d vector_cross (Vector3d & u, Vector3d & v) {

	Vector3d w ;
	w << 0, 0, 0 ;

	for (int i = 0 ; i < w.size() ; ++i) {

		for (int j = 0 ; j < w.size() ; ++j) {

			if (i == j) {
				// Do nothing.
			}

			else if (i == 0 && j == 1) {
				w(2) += u(i) * v(j) ;
			}

			else if (i == 0 && j == 2) {
				w(1) -= u(i) * v(j) ;
			}

			else if (i == 1 && j == 0) {
				w(2) -= u(i) * v(j) ;
			}

			else if (i == 1 && j == 2) {
				w(0) += u(i) * v(j) ;
			}

			else if (i == 2 && j == 0) {
				w(1) += u(i) * v(j) ;
			}

			else if (i == 2 && j == 1) {
				w(0) -= u(i) * v(j) ;
			}

		}

	}

	return w ;
}

double vector_dot (vector<double> & a, vector<double> & b) {

	if ( a.size() != b.size() ) {
		
		throw invalid_argument( "\n\n*****\nInput vector dimension problem in M2DO::FEA vector_dot ().\n*****\n\n" );
	
	}

	double dot = 0.0 ;

	for (int i = 0 ; i < a.size() ; ++i) {

		dot += a[i] * b[i] ;

	}

	return dot ;

}

VectorXd vector_to_eigen_vectorxd (vector<double> v) {

	Map<VectorXd> V (v.data(), v.size());
	
	return V ;

}


vector<double> eigen_vectorxd_to_vector(VectorXd V) {

	vector<double> v (V.data(), V.data() + V.rows() * V.cols()) ;
	
	return v ;

}

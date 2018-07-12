#include "linear_shape_function.h"

using namespace M2DO_FEA ;

LinearShapeFunction :: LinearShapeFunction () {

	spacedim = dim = -1 ;

}

LinearShapeFunction :: LinearShapeFunction (int spacedim, int dim) : spacedim (spacedim), dim (dim) {

	eta_values = MatrixXd::Constant(pow(2, spacedim), spacedim, -1) ;
	vector<int> eta_count (spacedim, 0) ;
	eta_count[0] += 1 ;

	for (int i = 1 ; i < pow(2, spacedim) ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			eta_values(i, j) = eta_values(i-1, j)  ;
			eta_count[j] += 1 ;

			if ( eta_count[j] == pow(2, max(1, j)) ) {
				eta_count[j] = 0 ;
				eta_values(i, j) *= -1  ;
			}

		}

	}

}

vector<double> LinearShapeFunction :: GetEta (int number) {

	/* 
		Returns the natural coordinates of desired shape function.
		
		This is quite a tricky loop and deserves more explanation.
		It assumes the "standard node numbering". Produces sequence like:
		-1 -1 -1
		+1 -1 -1
		+1 +1 -1
		-1 +1 -1
		-1 -1 +1
		... etc.
	*/

	vector<double> eta_vals  (spacedim, -1) ;
	vector<int>    eta_count (spacedim, 0) ;
	eta_count[0] += 1 ;

	for (int i = 0 ; i < number ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			eta_count[j] += 1 ;

			if ( eta_count[j] == pow(2, max(1, j)) ) {
				eta_count[j] = 0 ;
				eta_vals[j] *= -1 ;
			}

		}

	}

	return eta_vals ;

}

double LinearShapeFunction :: GetShapeFunctionGradients (int number, int component, vector<double> & eta) {

	double grad_val = 1.0 / pow(2, spacedim) ;

	for (int i = 0 ; i < spacedim ; ++i) {
		
		if (i == component) {
			grad_val *= eta_values(number, i) ;
		}
		
		else {
			grad_val *= 1 + eta_values(number, i) * eta[i] ;
		}

	}

	return grad_val ;

}

double LinearShapeFunction :: GetShapeFunctionValues (int number, vector<double> eta) {

	// For these linear elements, there are 2^spacedim shape functions.

	double val = 1.0 / pow(2, spacedim) ;
	vector<double> eta_vals = GetEta(number) ;

	for (int i = 0 ; i < spacedim ; ++i) {
		val *= 1 + eta_vals[i] * eta[i] ;
	}

	return val ;

}

VectorXd LinearShapeFunction :: GetShapeFunctionValuesVector (vector<double> eta) {

	// For these linear elements, there are 2^spacedim shape functions.

	VectorXd val_vec = VectorXd::Zero(pow(2, spacedim)) ;
	double val ;
	vector<double> eta_vals ;

	for (int number = 0 ; number < pow(2, spacedim); ++number) {

		val = 1.0 / pow(2, spacedim) ;
		eta_vals = GetEta(number) ;

		for (int i = 0 ; i < spacedim ; ++i) {
			val *= 1 + eta_vals[i] * eta[i] ;
		}

		val_vec(number) = val ;

	}

	return val_vec ;

}

VectorXd LinearShapeFunction :: GetShapeFunctionValuesFullVector (double v, int component) {

	/*

	Outputs a vector with extra zeros, for example, v_full = [v_0, 0].
	
	*/

	VectorXd v_full = VectorXd::Zero(dim) ;

	v_full(component) = v ;

	return v_full ;

}


VectorXd LinearShapeFunction :: GetShapeFunctionGradientsVector (int number, vector<double> & eta) {

	/*

	Same as for GetShapeFunctionGradients() above, though outputs all components as a vector.
	
	*/

	VectorXd grad_val_vec = VectorXd::Zero(spacedim) ;

	for (int k = 0 ; k < spacedim ; ++k) {

		double grad_val = 1.0 / pow(2, spacedim) ;

		for (int i = 0 ; i < spacedim ; ++i) {
			
			if (i == k) {
				grad_val *= eta_values(number, i) ;
			}
			
			else {
				grad_val *= 1 + eta_values(number, i) * eta[i] ;
			}

		}

		grad_val_vec(k) = grad_val ;

	}

	return grad_val_vec ;

}

VectorXd LinearShapeFunction :: GetShapeFunctionGradientsFullVector (VectorXd & v, int component) {

	/*

	Same as for GetShapeFunctionGradients_vec() above, though outputs a vector with extra zeros, for example,
	v_full = [grad(v_0), grad(0)].
	
	*/

	VectorXd v_full = VectorXd::Zero(spacedim * dim) ;

	int k_0 = spacedim * component ;
	v_full.segment(k_0, spacedim) = v ;

	return v_full ;

}



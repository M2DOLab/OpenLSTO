#include "quadrature.h"
#include "linear_shape_function.h"
#include "mesh.h"
#include "solid_element.h"

using namespace M2DO_FEA ;

SolidElement :: SolidElement (int spacedim, int order, Mesh & mesh) : spacedim (spacedim), order (order), mesh (mesh) {

	// Constructor fills node vector with -1 by default:
	node_ids = vector<int> (pow(2, mesh.spacedim), -1) ;
	material_id = 0 ;
	area_fraction = 1.0 ;
	linear_shape_function = LinearShapeFunction (spacedim, spacedim) ;
	quadrature            = GaussianQuadrature  (spacedim, order) ;

}

void SolidElement :: Print () {

	cout << "SolidElement (" ;

	for (int i = 0 ; i < node_ids.size() ; ++i) {

		if (i > 0) {

			cout << ", " ;

		}

		cout << node_ids[i] ;

	}

	cout << ")" ;

}

MatrixXd SolidElement :: J (vector<double> & eta) {

	MatrixXd J_mat = MatrixXd::Zero (spacedim, spacedim) ;

	for (int i = 0 ; i < spacedim ; ++i) {

		for (int j = 0 ; j < spacedim ; ++j) {

			for (int k = 0 ; k < pow(2, spacedim) ; ++k) {

				J_mat (i, j) += mesh.nodes[node_ids[k]].coordinates[j] * linear_shape_function.GetShapeFunctionGradients (k, i, eta) ;

			}

		}

	}

	return J_mat ;

}

MatrixXd SolidElement :: B (vector<double> & eta) {

	VectorXd shape_grad_j, shape_grad_j_full ;

	/*
		grad(u(x)) = [B] * {u}
		So we aim to find B here.
	*/

	MatrixXd B_mat = MatrixXd::Zero(spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int shape_j = 0, dim_j = 0 ;

	MatrixXd J_mat = J (eta) ;
	MatrixXd J_inv = J_mat.inverse() ;

	// Build the B matrix at the given eta point:
	for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

		shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
		shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;

		B_mat.col(j)      = shape_grad_j_full ;

		if (dim_j < spacedim-1) {
			dim_j = dim_j + 1 ;
		}

		else {
			dim_j = 0 ;
			shape_j = shape_j + 1 ;
		}

	} // for j (columns in B).

	return B_mat ;

}

MatrixXd SolidElement :: B_int () {

	MatrixXd J_mat, J_inv, B_int = MatrixXd::Zero (spacedim*spacedim, pow(2, spacedim) * spacedim) ;
	VectorXd shape_grad_j, shape_grad_j_full ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Add to the K matrix:
		*/

		B_int += B_mat * w * J_mat.determinant() ;

	} // for k (gauss points).

	return B_int ;

}

MatrixXd SolidElement :: K () {

	MatrixXd J_mat, J_inv, K_mat = MatrixXd::Zero (pow(2, spacedim) * spacedim, pow(2, spacedim) * spacedim) ;
	MatrixXd C = mesh.solid_materials[material_id].C ;
	VectorXd shape_grad_j, shape_grad_j_full ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Add to the K matrix:
		*/

		K_mat += B_mat.transpose() * C * B_mat * w * J_mat.determinant() ;

	} // for k (gauss points).

	return K_mat ;

}

MatrixXd SolidElement :: K_NodalProperties (VectorXd & nodal_properties, double min_property) {

	MatrixXd J_mat, J_inv, K_mat = MatrixXd::Zero (pow(2, spacedim) * spacedim, pow(2, spacedim) * spacedim) ;
	MatrixXd C1 = mesh.solid_materials[material_id].C ;
	VectorXd shape_grad_j, shape_grad_j_full ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Calculate the nodal property value at this gauss point:
		*/

		double nodal_properties_g = linear_shape_function.GetShapeFunctionValuesVector(eta).transpose() * nodal_properties ;

		/*
			Retrieve C from the material:
		*/

		MatrixXd C  = C1 + nodal_properties_g * (min_property*C1 - C1) ;

		/*
			Add to the K matrix:
		*/

		K_mat += B_mat.transpose() * C * B_mat * w * J_mat.determinant() ;

	} // for k (gauss points).

	return K_mat ;

}

vector<MatrixXd> SolidElement :: dKdz (double min_property) {

	MatrixXd J_mat, J_inv;
	vector<MatrixXd> dKdz_mat (node_ids.size(), MatrixXd::Zero (8, 8)) ;
	MatrixXd C1 = mesh.solid_materials[material_id].C ;
	VectorXd shape_grad_j, shape_grad_j_full ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Calculate shape function vector at this gauss point:
		*/

		VectorXd shape_vec = linear_shape_function.GetShapeFunctionValuesVector(eta);

		/*
			Add to the K matrix:
		*/

		for (int i = 0 ; i < node_ids.size() ; ++i) {

			MatrixXd dCdz  = shape_vec[i] * (min_property*C1 - C1) ;
			dKdz_mat[i] += B_mat.transpose() * dCdz * B_mat * w * J_mat.determinant() ;

		}

	} // for k (gauss points).

	return dKdz_mat ;

}

MatrixXd SolidElement :: M () {

	MatrixXd J_mat, J_inv, M_mat = MatrixXd::Zero (pow(2, spacedim) * spacedim, pow(2, spacedim) * spacedim) ;
	VectorXd shape_grad_j, shape_grad_j_full ;

	double rho = mesh.solid_materials[material_id].rho;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the H matrix at this gauss point: (2D only)
		*/

		MatrixXd H_mat = MatrixXd::Zero (2, 8) ;
		VectorXd     v = linear_shape_function.GetShapeFunctionValuesVector(eta) ;

		H_mat << v[0], 0.0, v[1], 0.0, v[2], 0.0, v[3], 0.0,
				  0.0, v[0], 0.0, v[1], 0.0, v[2], 0.0, v[3] ;

		// Update gauss point.
		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Add to the K matrix:
		*/

		M_mat += rho * H_mat.transpose() * H_mat * w * J_mat.determinant() ;

	} // for k (gauss points).

	return M_mat ;

}

MatrixXd SolidElement :: FThermalExpansion (double DeltaT) {

	MatrixXd J_mat, J_inv, f_thermal = MatrixXd::Zero (pow(2, spacedim) * spacedim, 1) ;
	MatrixXd C = mesh.solid_materials[material_id].C ;
	VectorXd shape_grad_j, shape_grad_j_full ;
	VectorXd epsilon = VectorXd::Zero(4);
	epsilon[0] = 1;
	epsilon[3] = 1;

	double alphaT = mesh.solid_materials[material_id].alphaT ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Add to the K matrix:
		*/

		f_thermal += B_mat.transpose() * C * epsilon * alphaT * DeltaT * w * J_mat.determinant() ;

	} // for k (gauss points).

	return f_thermal ;

}

MatrixXd SolidElement :: FThermalExpansion (VectorXd deltaT_vector) {

	MatrixXd J_mat, J_inv, f_thermal = MatrixXd::Zero (pow(2, spacedim) * spacedim, 1) ;
	MatrixXd C = mesh.solid_materials[material_id].C ;
	VectorXd shape_grad_j, shape_grad_j_full ;
	VectorXd epsilon = VectorXd::Zero(4);
	epsilon[0] = 1;
	epsilon[3] = 1;

	double alphaT = mesh.solid_materials[material_id].alphaT ;

	double w ;
	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	/*
		grad(u(x)) = [B] * {u}
	*/

	MatrixXd B_mat = MatrixXd::Zero (spacedim * spacedim, pow(2, spacedim) * spacedim) ;

	int n_gauss = pow (order, spacedim) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		int shape_j = 0, dim_j = 0 ;

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc. Same goes
			for the weighting w.
		*/

		w = 1.0 ;

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;
			w      *= quadrature.w[eta_count[l]] ;

		}


		VectorXd lsfvals =  linear_shape_function.GetShapeFunctionValuesVector(eta);

		double DeltaT = deltaT_vector.dot(lsfvals);

		J_mat = J (eta) ;
		J_inv = J_mat.inverse() ;

		/*
			Build the B matrix at this gauss point:
		*/

		for (int j = 0 ; j < spacedim * pow(2, spacedim) ; ++j) {

			shape_grad_j      = J_inv * linear_shape_function.GetShapeFunctionGradientsVector(shape_j, eta) ;
			shape_grad_j_full = linear_shape_function.GetShapeFunctionGradientsFullVector(shape_grad_j, dim_j) ;


			B_mat.col(j)      = shape_grad_j_full ;

			if (dim_j < spacedim-1) {
				dim_j = dim_j + 1 ;
			}

			else {
				dim_j = 0 ;
				shape_j = shape_j + 1 ;
			}

		} // for j (columns in B).

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

		/*
			Add to the K matrix:
		*/

		f_thermal += B_mat.transpose() * C * epsilon * alphaT * DeltaT * w * J_mat.determinant() ;

	} // for k (gauss points).

	return f_thermal ;

}

VectorXd SolidElement :: NaturalToPhysicalCoordinates (vector<double> & eta) {

	VectorXd x = VectorXd::Zero (spacedim) ;

	VectorXd v = linear_shape_function.GetShapeFunctionValuesVector (eta) ;

	for (int i = 0 ; i < spacedim ; ++i) {

		for (int j = 0 ; j < v.size() ; ++j) {

			x[i] += v[j] * mesh.nodes[node_ids[j]].coordinates[i] ;

		}

	}

	return x ;

}

MatrixXd SolidElement :: PhysicalGaussPoissCoordinates () {

	int n_gauss = pow (order, spacedim) ;

	MatrixXd X = MatrixXd::Zero (n_gauss, spacedim) ;

	vector<double> eta (spacedim, 0), eta_count (spacedim, 0) ;

	for (int k = 0 ; k < n_gauss ; ++k) {

		/*
			The first gauss point is located at eta = [quadrature.eta[0], quadrature.eta[0], ...].
			Since eta_count was initialized with zeros in all dimensions, we don't need to do
			anything fancy yet; just set eta[l] = quadrature.eta[eta_count[l]], etc.
		*/

		for (int l = 0 ; l < spacedim ; ++l) {

			eta[l]  = quadrature.eta[eta_count[l]] ;

		}

		X.row(k) = NaturalToPhysicalCoordinates(eta) ;

		eta_count = quadrature.UpdateEtaCounter(eta_count) ;

	} // for k (gauss points).

	return X ;

}

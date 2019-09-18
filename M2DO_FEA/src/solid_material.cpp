#include "solid_material.h"

using namespace M2DO_FEA ;

SolidMaterial :: SolidMaterial (int spacedim, double E, double nu, double rho, double alphaT, double h) : spacedim (spacedim), E (E), nu (nu), rho (rho), alphaT (alphaT), h (h) {

	if (spacedim == 2) {

		MatrixXd A = MatrixXd::Zero (4, 4) ;

		A << 1,  0,   0,  0,
			 0, 0.5, 0.5, 0,
			 0, 0.5, 0.5, 0,
			 0,  0,   0,  1 ;

		// Voight matrix:
		V = MatrixXd::Zero (4, 4) ;

		V <<    1,   0,   0, -0.5,
			    0, 1.5,   0,    0,
			    0,   0, 1.5,    0,
			 -0.5,   0,   0,    1 ;

		// Note: This is the plane stress formulation!
		MatrixXd D = MatrixXd::Zero (4, 4) ;

		D << 1,     0, 		  0,    nu,
			 0, (1-nu)/2, (1-nu)/2, 0,
			 0, (1-nu)/2, (1-nu)/2, 0,
			 nu,     0,       0,    1 ;

		D *= E / (1-pow(nu,2)) ;

		C = h * D * A ;

		// Note: This is the axisymmetric formulation!
		MatrixXd D_axisymmetric = MatrixXd::Zero (5, 5) ;
		
		D_axisymmetric << 1-nu,     0, 		 	0,           nu, 	nu,
							 0, 	(1-2*nu)/2, (1-2*nu)/2,  0, 	 0,
							 0, 	(1-2*nu)/2, (1-2*nu)/2,  0,  	 0,
						    nu,		0,       	0,  		 1-nu,   nu,
					        nu, 	0,       	0,  		 nu,     1-nu;
		
		D_axisymmetric *= E / ((1 + nu)*(1-2*nu)) ;

		C_axisymmetric = h * D_axisymmetric ;

	}

	else if (spacedim == 3) {

		MatrixXd A = MatrixXd::Zero (9, 9) ;

		A << 1,   0,   0,   0, 0,   0,   0,  0,  0,
			 0, 0.5,   0, 0.5, 0,   0,   0,  0,  0,
			 0,   0, 0.5,   0, 0,   0, 0.5,  0,  0,
			 0, 0.5,   0, 0.5, 0,   0,   0,  0,  0,
			 0,   0,   0,   0, 1,   0,   0,  0,  0,
			 0,   0,   0,   0, 0, 0.5,   0, 0.5, 0,
			 0,   0, 0.5,   0, 0,   0, 0.5,  0,  0,
			 0,   0,   0,   0, 0, 0.5,   0, 0.5, 0,
			 0,   0,   0,   0, 0,   0,   0,  0,  1 ;

		MatrixXd D = MatrixXd::Zero (9, 9) ;

		D << (1-nu),        0,        0,        0,     nu,        0,        0,        0,      nu,
			      0, (1-2*nu),        0,        0,      0,        0,        0,        0,       0,
			      0,        0, (1-2*nu), 	    0,      0,        0,        0,        0,       0,
			      0,	    0, 		  0, (1-2*nu),      0,        0,        0,        0,       0,
			     nu, 	    0,        0,        0, (1-nu),        0,        0,        0,       nu,
			      0,        0,        0,        0,      0, (1-2*nu),        0,        0,       0,
				  0,        0,        0,        0,      0,        0, (1-2*nu),        0,       0,
			      0,        0,        0,        0,      0,        0,        0, (1-2*nu),       0,
			     nu,        0,        0,        0,      nu,        0,        0,        0,  (1-nu) ;

		D *= E / ((1+nu) * (1-2*nu)) ;

		C = D * A ;

		// Voight matrix:
		V = Eigen::Matrix<double, 9, 9>::Zero() ;
	 V << 1,   0,   0,   0, -0.5,   0,   0,  0,  -0.5,
			0, 1.5,   0, 0.0, 0,   0,   0,  0,  0,
			0,   0, 1.5,   0, 0,   0, 0.0,  0,  0,
			0, 0.0,   0, 1.5, 0,   0,   0,  0,  0,
			-0.5,   0,   0,   0, 1,   0,   0,  0,  -0.5,
			0,   0,   0,   0, 0, 1.5,   0, 0.0, 0,
			0,   0, 0.0,   0, 0,   0, 1.5,  0,  0,
			0,   0,   0,   0, 0, 0.0,   0, 1.5, 0,
			-0.5,   0,   0,   0, -0.5,   0,   0,  0,  1 ;

	}

}

void SolidMaterial :: Print () {

	std::cout << "Solid Material (E = " << E << ", nu = " << nu << ", rho = " << rho << ", h = " << h << ")" ;

}

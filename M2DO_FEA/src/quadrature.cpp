#include "quadrature.h"

using namespace M2DO_FEA ;

GaussianQuadrature :: GaussianQuadrature () {

	spacedim = order = -1 ;

}

GaussianQuadrature :: GaussianQuadrature (int spacedim, int order) : spacedim (spacedim), order (order), eta (order, 0), w (order, 0) {

	if (order == 1) {

		eta[0] = 0 ;

		w[0] = 2 ;

	}

	else if (order == 2) {

		eta[0] = -1.0 / sqrt(3) ;
		eta[1] = 1.0 / sqrt(3) ;

		w[0] = 1.0 ;
		w[1] = 1.0 ;

	}

	else if (order == 3) {

		eta[0] = -sqrt(3.0 / 5) ;
		eta[1] = 0 ;
		eta[2] = sqrt(3.0 / 5) ;

		w[0] = 5.0 / 9 ;
		w[1] = 8.0 / 9 ;
		w[2] = 5.0 / 9 ;
		
	}

	else if (order == 4) {

		eta[0] = -sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)) ;
		eta[1] = -sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)) ;
		eta[2] =  sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)) ;
		eta[3] =  sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)) ;

		w[0] = (18.0 - sqrt(30.0))/36 ;
		w[1] = (18.0 + sqrt(30.0))/36 ;
		w[2] = (18.0 + sqrt(30.0))/36 ;
		w[3] = (18.0 - sqrt(30.0))/36 ;
		
	}

	else if (order == 5) {

		eta[0] = -1.0/3 * sqrt(5.0 + 2.0*sqrt(10.0/7)) ;
		eta[1] = -1.0/3 * sqrt(5.0 - 2.0*sqrt(10.0/7)) ;
		eta[2] =  0.0 ;
		eta[3] =  1.0/3 * sqrt(5.0 - 2.0*sqrt(10.0/7)) ;
		eta[4] =  1.0/3 * sqrt(5.0 + 2.0*sqrt(10.0/7)) ;

		w[0] = (322.0 - 13.0*sqrt(70.0)) / 900 ;
		w[1] = (322.0 + 13.0*sqrt(70.0)) / 900 ;
		w[2] = 128.0 / 155 ;
		w[3] = (322.0 + 13.0*sqrt(70.0)) / 900 ;
		w[4] = (322.0 - 13.0*sqrt(70.0)) / 900 ;
		
	}

}


void GaussianQuadrature :: Print () {

	cout << "GaussianQuadrature( eta(" ;

	for (int i = 0 ; i < order ; ++i) {

		if (i > 0) {
			cout << ", " ;
		}

		cout << eta[i] ;
	}

	cout << "), w(" ;

	for (int i = 0 ; i < order ; ++i) {

		if (i > 0) {
			cout << ", " ;
		}

		cout << w[i] ;
	}

	cout << ") )" ;
}


vector<double> GaussianQuadrature :: UpdateEtaCounter (vector<double> & eta_count) {

	/*
	We now update eta_count, which is used in the next iteration to find
	the eta coordinates of the next gauss point. It is actually a tricky
	loop, because for example if spacedim = 3, you want to produce 
	a sequence of eta_count like:

	0	0	0
	1	0	0
	2	0	0
	0	1	0
	1	1	0
	2	1	0
	0	2	0
	1	2	0
	2	2	0
	0	0	1
	... etc.

	The following code will produce such a sequence.
	*/

	eta_count[0] += 1 ;

	if (eta_count[0] > (order - 1)) {

		eta_count[0] = 0 ;
		
		for (int l = 1 ; l < spacedim ; ++l) {

			eta_count[l] += 1 ;

			if (eta_count[l] <= (order - 1)) {
				break ;
			}

			else {
				eta_count[l] = 0 ;
			}

		}

	}

	return eta_count ;
}


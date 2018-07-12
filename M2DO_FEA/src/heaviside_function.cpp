#include "heaviside_function.h" 

using namespace M2DO_FEA ;

HeavisideFunction :: HeavisideFunction () : delta (0.0), beta (0.0) {

	//
	
}

HeavisideFunction :: HeavisideFunction (double delta, double beta) : delta (delta), beta (beta) {

	//
	
}

void HeavisideFunction :: print () {

	cout << "HeavisideFunction (delta = " << delta << ", beta = " << beta << ")" ;

}

double HeavisideFunction :: value (double x) {

	x -= beta ;

	if (x <= -delta) {

		return 1.0 ;

	}

	else if (x >= delta) {

		return 0.0 ;

	}

	else {

		double z = (x+delta)/(2.0*delta) ;
		double y = 1.0 - 6.0*pow(z, 5.0) + 15.0*pow(z, 4.0) - 10.0*pow(z, 3.0) ;
		// double y = 1.0 + 20.0*pow(z, 7.0) - 70.0*pow(z, 6.0) + 84.0*pow(z, 5.0) - 35.0*pow(z, 4.0) ;
		// double y = 1.0 + 3432.0*pow(z, 15.0) - 25740.0*pow(z, 14.0) + 83160.0*pow(z, 13.0) - 150150.0*pow(z, 12.0) + 163800.0*pow(z, 11.0) - 108108.0*pow(z, 10.0) + 40040.0*pow(z, 9.0) - 6435.0*pow(z, 8.0) ;
		return y ;

	}

}

double HeavisideFunction :: grad (double x) {

	x -= beta ;

	if (x <= -delta) {

		return 0.0 ;

	}

	else if (x >= delta) {

		return 0.0 ;

	}

	else {

		double    z = (x+delta)/(2.0*delta) ;
		double dzdx = 1.0/(2.0*delta) ;
		double dydz = -30.0*pow(z, 4.0) + 60.0*pow(z, 3.0) - 30.0*pow(z, 2.0) ;
		// double dydz = 140.0*pow(z, 6.0) - 420.0*pow(z, 5.0) + 420.0*pow(z, 4.0) - 140.0*pow(z, 3.0) ;
		// double dydz = 51480.0*pow(z, 14.0) - 360360.0*pow(z, 13.0) + 1081080.0*pow(z, 12.0) - 1801800.0*pow(z, 11.0) + 1801800.0*pow(z, 10.0) -  1081080.0*pow(z, 9.0) + 360360.0*pow(z, 8.0) - 51480.0*pow(z, 7.0) ;
		return (dydz * dzdx) ;

	}

}
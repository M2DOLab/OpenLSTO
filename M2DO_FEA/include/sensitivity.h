#ifndef M2DO_FEA_SENSITIVITY_H
#define M2DO_FEA_SENSITIVITY_H

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

#include <../../vendor/eigen3/Eigen/Dense>
#include <../../vendor/eigen3/Eigen/Sparse>

using namespace std ;
using namespace Eigen ;

namespace M2DO_FEA {

    // -------------------------------------------------------------------------------- //
    // LEAST SQUARES EXTERN FUNCTION
    extern "C" int dgels_(const char *trans, int *m, int *n, int *nrhs,
            double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

    // -------------------------------------------------------------------------------- //
    // SENSITIVITIES CLASS
    class Sensitivity {

        private:
            //

        public:
            // Sensitivities at prescribed point (is a vector so we can compute at Gauss points when handy)
            vector<double> sensitivity_at_gauss_point;
	    vector<double> sensitivity_component1_at_gauss_point;
	    vector<double> sensitivity_component2_at_gauss_point;
	    vector<double> sensitivity_component3_at_gauss_point;
            // Global coordinates of the sensitivities points
            vector< std::vector<double> > coordinate;

	    // Average sensivity (calculated from gauss point sensitivities)
	    double sensitivity_average;

            // Is element excluded from the optimisation?
           bool isExcluded;


    };

    // STRESS, STRAIN AND ADJOINT STRAIN CLASS
    class StressStrain {

        private:
            //

        public:
            // von Mises stress.
            std::vector<double> von_mises; // For all Gauss points.
            double von_mises_average; // Average for the element.

            // Strain.
            std::vector<double> strain; // For all Gauss points.

    };


    // -------------------------------------------------------------------------------- //
    // LEAST SQUARES CLASS
    class LeastSquares {
        public:
            // Distances from gauss point to boundary point
            double distance_from_gauss_point;
            // Area fraction at the Gauss point
            double area_fraction_at_gauss_point;
            // Element number
            int element_number;
            // Gauss point local number
            int gauss_point_number;
            // Coordinates of the Gauss point
            vector<double> coordinate;

    };

    // -------------------------------------------------------------------------------- //
    // SENSITIVITY Analysis CLASS
    class SensitivityAnalysis {

        private:
            // Compute Gauss points global coordinates.
            void ComputeSensitivitiesCoordinates (bool time_it) ;

	    // Solve least squares problem.
            double SolveLeastSquares(vector<LeastSquares> least_squares, vector<double> boundary_point, int indicator = 0);

            std::vector<double> mat_vec_mult( std::vector<std::vector<double>> &AtA, std::vector<double> &v_in );

            double vec_vec_mult( std::vector<double> &v_in1, std::vector<double> &v_in2 );

            // Least squares information class.
            vector<LeastSquares> least_squares;

        public:
            // Properties:

            int spacedim, order, dim;

            // Vector of sensitivities.
            vector<Sensitivity> sensitivities;
            // Boundary Sensitivities.
            vector<double> boundary_sensitivities;

            // Attaching classes.
            StationaryStudy & study;

            // Vector of strains.
            vector<StressStrain> strains;

            // Maximum von Mises stress.
            double von_mises_max;


            // For stress analysis.
            double objective; // Objective function value.

            // General methods:
            SensitivityAnalysis(StationaryStudy & study); // Constructor

            // Exclude elements from sensitivity computation.
            void ExcludeElements (vector<int> selected_elements);

            // Compliance problem.
            void ComputeComplianceSensitivities (bool time_it);


            // Least squares functionalities.
            void ComputeBoundarySensitivities(vector<double> boundary_point, double radius = 2, int indicator = 0, double p_norm = 6);

            // Stress problem
            void ComputeStressSensitivities (bool time_it, double pnorm);

            // stress problem for 3D
            void ComputeStressSensitivities3D (bool time_it, double pnorm);





	    // Printing
	    void WriteAverageVonMisesTxt (int datapoint, int num_elem_x, int num_elem_y, std::string filePath = "", std::string txtFileName = "stress");
	    void WriteAverageVonMisesVtk ();
	    void WriteAverageSensitivitiesTxt (int datapoint, int num_elem_x, int num_elem_y, std::string filePath = "", std::string txtFileName = "sens");
    };

}

#endif

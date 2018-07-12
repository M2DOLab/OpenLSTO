#ifndef LSM_OPTI_3D_H
#define LSM_OPTI_3D_H

#include <cmath>
#include <vector>



namespace FEA = M2DO_FEA ;

typedef unsigned int uint;

using namespace std;
using namespace Eigen ;

class SensitivityData
{
  /*
    This sensitivity data class can be used to pass info between FE and LSM meshes
  */
  public:

  double compliance = 0; // compliance
  double Vol; // Volume of solid
  int iter; // iteration number
  int algo = 0; // algo = 0 for newton-raphson
  int nx; // num of lsm elements in x direction
  int ny; // num of lsm elements in y direction
  int nz; // num of lsm elements in z direction

  int nx1, nz1; //  parameters of the L beam

  std::vector<double> volumeFractions; // vector of volume fraction
  std::vector<double> bsens; // boundary sensitivity
  std::vector<double> opt_vel_simplex; // optimum velocities computed using simplex
  std::vector<double> opt_vel_nlopt; // optimum velocities computed using nlopt
  std::vector<double> pointAreas;// boundary point areas
  std::vector<double> bPoints; // boundary points
  std::vector<double> uguess; // guess solution

  std::vector<double> vsens; // volume sensitivity

  // bounda for using simplex
  std::vector<double> LB;
  std::vector<double> UB;
  std::vector<double> LAM;

  // maps ls and fe meshes
  std::vector<double> LS2FEmap;


  double move_limit; // cfl condition



  double MaxVol; // maximum volume in percent

  int bpointsize; // number of boundary points

  double sensi_time;// time taken by sensitivity
  double solver_time;// time taken by sensitivity

};


void PerformOptimization (SensitivityData &SensData); // optimizes the boundary velocities

void PerformOptimization_LBeam (SensitivityData &SensData); // optimized the boundary velocities for the L beam

#include "lsm_opti_3d.cpp"

#endif

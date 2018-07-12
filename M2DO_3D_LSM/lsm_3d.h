#ifndef LSM_3D_H
#define LSM_3D_H


#include <cmath>
#include <vector>
#include <algorithm>
#include "marching_cubes_cross.h"
#include "mp_vector.h"

#include "marching_cubes_cross.cpp"
#include "mp_vector.cpp"



typedef unsigned int uint;

using namespace std;

class LevelSet3D
{
  /*
    This is the first version of the M2DO level set module.
    This class members are:
      the 3d grid dimensions: nx, ny, nz
      signed distance on the grid points: an mpVector phi
      status of a point (whether it is inside or outside of a narrow band)
      boundary points
      boundary triangles

  */
  public:

  uint num_elem_x, num_elem_y, num_elem_z; // Number of elements in x,y,z directions

  uint num_grid_pts; // Number of grid points

  uint num_cells; // Number of cells

  std::vector< std::vector<double> > boundary_pts; // Vector of boundary pts

  std::vector<double>  boundary_pts_one_vector; // Vector of boundary pts in a single vector

  std::vector<double> volumefraction_vector ; // Vector of volume fractions in the xyz axis

  std::vector<double> volumefraction_vector_perturb ; // Vector of volume fractions used for pertrubations

  std::vector<double> boundary_areas; // Vector of boundary areas

  std::vector<double> opt_vel; // vector of optimum velocities

  uint num_boundary_pts; // Number of boundary points

  mp4Vector *phi ; // pointer to signed distance mp4vector

  std::vector<double> grid_vel; // grid velocities which are extrapolated from the optimum velocities

  std::vector<double> grid_gradient; // grid gradients

  float iso_value = 0.0;    // value at which surface occurs

  TRIANGLE *triangle_array; // pointer to triangle array

  int num_triangles;        // number of triangles in triangle array

  std::vector<int> indices_considered_inside; //  indices of values considered inside

  std::vector<int> indices_considered_outside; //  indices of values considered inside

  std::vector<int> indices_considered; //  indices of points considered in FastMarchingMethod

  std::vector<double> phi_considered; //  phi of points considered in FastMarchingMethod

  std::vector<int> phi_status; //  status of a point: whether it is considered or not

  std::vector<int> vol_frac_status; //  status of a point: whether it is considered or not

  std::vector<double> phi_temp; // a temperory phi vector used in setting up the narrow band

  int narrow_band_width = 3;

  int num_gauss_pts = 2;

  std::vector<std::vector<double> > holes; // holes

  std::vector<std::vector<double> > cubic_holes; // holes


  /*
    Functions
  */
  void SetBoxDimensions(uint box_x, uint box_y, uint box_z); // Sets box dimensions

  uint nx1, nz1;

  void MakeBox(); // Initializes the signed distance function for a box;

  void MakeLBeam(); // Initializes an L beam

  void MarchingCubesWrapper(); //  returns faces, midpoints, and vertices of the discretized boundary

  void SetupNarrowBand(); // sets up a narrow band around the boundary points

  void ExtrapolateVelocities(); //  extrapolates velocities from boundary points to grid points

  void FastMarchingMethod(); // solves the Eikonal equation to compute phi

  int Grid_pt_to_index_zyx(int x, int y, int z); // returns the index of a grid point in the zyx system (mp4vector convention)

  std::vector<int>  Index_to_grid_pt_zyx(int index); // returns the grid point of an index the zyx system (mp4vector convention)

  void Advect(); // advects the level set

  void Advect_LBeam(); // advects the level set for the Lbeam problem

  void ComputeGradients(); // computes gradients of grid points in the narrow band

  double gradHJWENO(double v1, double v2, double v3, double v4, double v5) ; // computes gradient

  void CalculateVolumeFractions(); // calculates the volume fractions in the xyz axis system

  void SolveEikonal(std::vector<int> indices_xyz); // solves the Eikonal equation for the point index_x, index_y, index_z

  void UpdateVelocity(int xin, int yin, int zin);// updates velocity

  void WriteFVfile();

  // write to stl
  void WriteSTL(int box_smooth);
};



#include "lsm_3d.cpp"

#endif

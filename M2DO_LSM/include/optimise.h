
  /*
    A class for finding the solution for the optimum velocity vector.
  */

#ifndef _OPTIMISE_H
#define _OPTIMISE_H


class Optimise
{
public:
    // Constructor.
    Optimise(std::vector<BoundaryPoint>&, double&, double &);


    double boundary_area; // area of structure
    double mesh_area; // area of the entire mesh
    double max_area; // maximum area
    double length_x; // length of structure in x direction
    double length_y; // length of structure in y direction

    void Solve_With_NewtonRaphson(); // function to calculate lambda_f

    void get_lambdas(std::vector<double> &lambdas);

private:
    /// The number of boundary points.
    unsigned int nPoints;

    /// A reference to a vector of boundary points.
    std::vector<BoundaryPoint>& boundaryPoints;

    /// The effective time step.
    double& timeStep;

    double & moveLimit; // moveLimit (CFL condition)

    double lambda_f;// psuedo lagrange parameter corresponding to objective
    double lambda_g;// psuedo lagrange parameter corresponding to constraint


};

#include "../src/optimise.cpp"

#endif  /* _OPTIMISE_H */

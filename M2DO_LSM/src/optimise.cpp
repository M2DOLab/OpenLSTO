

  /*
    A class for finding the solution for the optimum velocity vector.
  */


Optimise::Optimise(std::vector<BoundaryPoint>& boundaryPoints_,
                   double& timeStep_,
                   double& moveLimit_) :
                   boundaryPoints(boundaryPoints_),
                   timeStep(timeStep_),
                   moveLimit(moveLimit_)
{
  // Constructor
}


void Optimise::Solve_With_NewtonRaphson ()
{

    /*
      This method computes the optimum velocities of the boundary points using the Newton Raphson method

      We try to solve for lambda_f, such that

      bpointvelocity = lambda_f*sf + lambda_g*sg

      sf and sg are sensitivities corresponding to compliance and volume constraint
    */

      nPoints = boundaryPoints.size(); // define number of points

      timeStep = 1.0; // timeStep ,  DON'T CHANGE THIS

      int bpointsize = nPoints; // same as number of points

      double abssens = std::abs(boundaryPoints[0].sensitivities[0]); // absolute value of sensitivitiy initialized

      // normalize sensitivities of compliance
      for (int i = 0 ; i < bpointsize ; ++i) abssens = std::max(abssens, std::abs(boundaryPoints[i].sensitivities[0]));

      for (int i = 0 ; i < bpointsize ; ++i) boundaryPoints[i].sensitivities[0] == boundaryPoints[i].sensitivities[0]/abssens;


      // Set up the Newton-Raphson iteration problem
      // We solve this equation for lamda_f using Newton-Raphson
      // target_area = boundary_area + sum( cg_i(lambda_f*sf_i + lambda_g*sg_i ) )
      // with lambda_g = movelimit
      // target_area is the area of the structure that the optimizer will be aiming to target during this iteration

      lambda_g = moveLimit;

      double fraction_area = 0.5; // delta_area will be 50 percent of maximum area
      double target_area = boundary_area;
      for(int i = 0; i<bpointsize; i++){
        target_area += boundaryPoints[i].length*fraction_area*(-lambda_g);
      }
      target_area = std::max(max_area*mesh_area, target_area);

      std::vector<double> OptVel(bpointsize); // Initialise optimum velocity vector
      std::vector<double> curpt = {0,0};

      // compute domain distance
      // domain distance is the distace from the currt point to the domain distance
      std::vector<double> domain_distance_vector(bpointsize,0.0);
      for (int i = 0; i < bpointsize; i++)
      {
        curpt[0] = boundaryPoints[i].coord.x;
        curpt[1] = boundaryPoints[i].coord.y;

        double domdist; // distance from fem domain

        domdist = std::min({ std::abs( curpt[0] - 0.0),std::abs( curpt[0] - length_x),std::abs( curpt[1] - 0.0),std::abs( curpt[1] - length_y) });

        if(  ( curpt[0] - length_x >=0.0 || -(curpt[0] -0) >=0.0  || curpt[1] - length_y >=0.0 || -(curpt[1] -0) >=0.0  ) )
        {
          domdist = -1.0*domdist;
        }
        domain_distance_vector[i] = domdist;

      }

      double lambda_0 = 0.0; // initial guess for lambda_f
      double delta_lambda = 0.1; // delta_lambda to compute derivatives
      double lambda_cur;
      double new_area;

      double default_area = boundary_area;
      for(int i = 0; i<bpointsize; i++){
        default_area += boundaryPoints[i].length*std::min( domain_distance_vector[i], lambda_g*boundaryPoints[i].sensitivities[1] + lambda_0*boundaryPoints[i].sensitivities[0] );
      }

      // Newton-Raphson iterations start here
      for (int iter_NR = 0; iter_NR < 50; iter_NR ++)
      {
        lambda_cur = lambda_0 + 0*delta_lambda;
        // compute new area
        double new_area0 = boundary_area;
        for(int i = 0; i<bpointsize; i++){
          new_area0 += boundaryPoints[i].length*std::min( domain_distance_vector[i], lambda_g*boundaryPoints[i].sensitivities[1] + lambda_cur*boundaryPoints[i].sensitivities[0] );
        }

        lambda_cur = lambda_0 + delta_lambda;
        // compute new area
        double new_area2 = boundary_area;
        for(int i = 0; i<bpointsize; i++){
          new_area2 += boundaryPoints[i].length*std::min( domain_distance_vector[i], lambda_g*boundaryPoints[i].sensitivities[1] + lambda_cur*boundaryPoints[i].sensitivities[0] );
        }

        lambda_cur = lambda_0 - delta_lambda;
        // compute new area
        double new_area1 = boundary_area;
        for(int i = 0; i<bpointsize; i++){
          new_area1 += boundaryPoints[i].length*std::min( domain_distance_vector[i], lambda_g*boundaryPoints[i].sensitivities[1] + lambda_cur*boundaryPoints[i].sensitivities[0] );
        }

        double slope = (new_area2 - new_area1)/ 2 / delta_lambda;

        lambda_0 -= (new_area0 - target_area)/slope; // Newton-Raphson iterative equation

        if(std::abs(new_area0 - target_area) < 1.0E-3) break;
      }

      lambda_f = lambda_0;

      // assign the optimum velocities
      for (int i = 0; i < bpointsize; i++)
      {
        double domdist; // distance from fem domain

        domdist = domain_distance_vector[i];

        boundaryPoints[i].velocity = -1.0*std::min( lambda_f*boundaryPoints[i].sensitivities[0] + lambda_g*boundaryPoints[i].sensitivities[1], domdist);

      }

      // cfl condition: if abs(velocity) > moveLimit then scale velocities appropriately
      double absvel = 0.0;
      for (int i = 0; i < bpointsize; i++) absvel = std::max(absvel , boundaryPoints[i].velocity);

      if (absvel > moveLimit) for (int i = 0; i < bpointsize; i++) boundaryPoints[i].velocity = moveLimit*boundaryPoints[i].velocity/absvel;


}

void Optimise::get_lambdas(std::vector<double> &lambdas) {

  lambdas.clear();
  lambdas.push_back(-lambda_f);
  lambdas.push_back(-lambda_g);

}

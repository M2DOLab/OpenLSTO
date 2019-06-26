
void LevelSet3D::CalculateVolumeFractions()
{
  /*
    This funtion calculates the volume fraction vector in the xyz axis system
  */

  volumefraction_vector.resize(num_cells);


  #pragma omp parallel
  {
  #pragma omp for nowait
  for(int k = 0; k < num_elem_z; k++)
  {
    for(int j = 0; j < num_elem_y; j++)
    {
      for(int i = 0; i < num_elem_x; i++)
      {
        // The coordinates of the point in the lower left aft are (i,j,k)

        int counter_cell = i + num_elem_x*j + num_elem_x*num_elem_y*k ;

        // Save the indces of all the corner points
        std::vector<int> indices_temp(8,1);
        indices_temp[0] = Grid_pt_to_index_zyx(i,j,k);
        indices_temp[1] = Grid_pt_to_index_zyx(i+1,j,k);

        indices_temp[2] = Grid_pt_to_index_zyx(i,j+1,k);
        indices_temp[3] = Grid_pt_to_index_zyx(i+1,j+1,k);

        indices_temp[4] = Grid_pt_to_index_zyx(i,j,k+1);
        indices_temp[5] = Grid_pt_to_index_zyx(i+1,j,k+1);

        indices_temp[6] = Grid_pt_to_index_zyx(i,j+1,k+1);
        indices_temp[7] = Grid_pt_to_index_zyx(i+1,j+1,k+1);

        // Check if all the the phi values at the indices have the same sign
        // If not, then the element is cut --> volume fraction is in between 0 and 1

        if(phi[indices_temp[0]].val >= 0 && phi[indices_temp[1]].val >= 0 && phi[indices_temp[2]].val >= 0 && phi[indices_temp[3]].val >= 0 && phi[indices_temp[4]].val >= 0&& phi[indices_temp[5]].val >= 0 && phi[indices_temp[6]].val >= 0 && phi[indices_temp[7]].val >= 0)
        {
          // cell is completely inside
          volumefraction_vector[counter_cell] = 1.0;
        }
        else
        {
          if(phi[indices_temp[0]].val < 0 && phi[indices_temp[1]].val < 0 && phi[indices_temp[2]].val < 0 && phi[indices_temp[3]].val < 0 && phi[indices_temp[4]].val < 0&& phi[indices_temp[5]].val < 0 && phi[indices_temp[6]].val < 0 && phi[indices_temp[7]].val < 0)
          {
            // cell is completely outside
            volumefraction_vector[counter_cell] = 0.0;
          }
          else
          {
            // cell is cut
            // use gauss points to compute volume

            std::vector<double> weights(num_gauss_pts,1.0/(1.0*num_gauss_pts));
            std::vector<double> gaussPoints(num_gauss_pts,0.0);

            for(int i1 = 1; i1<= num_gauss_pts; i1++ ) gaussPoints[i1-1] = 1.0*i1/(1.0 + num_gauss_pts);

            double volume_fraction_temp = 0.0; // set volume fraction to zero

            // Loop through all gauss points
            for (uint gauss_i = 0; gauss_i < num_gauss_pts; gauss_i++) {
              for (uint gauss_j = 0; gauss_j < num_gauss_pts; gauss_j++) {
                for (uint gauss_k = 0; gauss_k < num_gauss_pts; gauss_k++) {

                  // Coordinates of gauss point in element space (local)
                  double xLocal = gaussPoints[gauss_i];
                  double yLocal = gaussPoints[gauss_j];
                  double zLocal = gaussPoints[gauss_k];

                  // Coordinate of gauss point in grid space (global)
                  double xGlobal = i + xLocal;
                  double yGlobal = j + yLocal;
                  double zGlobal = k + zLocal;

                  double N1 = phi[indices_temp[0]].val*(1-xLocal)*(1-yLocal)*(1-zLocal);
                  double N2 = phi[indices_temp[1]].val*(xLocal)*(1-yLocal)*(1-zLocal);

                  double N3 = phi[indices_temp[2]].val*(1-xLocal)*(yLocal)*(1-zLocal);
                  double N4 = phi[indices_temp[3]].val*(xLocal)*(yLocal)*(1-zLocal);

                  double N5 = phi[indices_temp[4]].val*(1-xLocal)*(1-yLocal)*(zLocal);
                  double N6 = phi[indices_temp[5]].val*(xLocal)*(1-yLocal)*(zLocal);

                  double N7 = phi[indices_temp[6]].val*(1-xLocal)*(yLocal)*(zLocal);
                  double N8 = phi[indices_temp[7]].val*(xLocal)*(yLocal)*(zLocal);

                  //double up = up1+up2;

                  double up = N1 + N2 + N3 + N4+ N5+ N6 + N7 + N8;

                  double u = 0.0;

                  if (up > 0.0 ) u = 1.0;

                  if(up <= 1.0/(num_gauss_pts) && up >= 0.0) u = (1.0*num_gauss_pts+1.0)*(up);

                  volume_fraction_temp += weights[gauss_i]*weights[gauss_j]*weights[gauss_k]*u;



                }
              }
            }// end of gauss points

            volumefraction_vector[counter_cell] = volume_fraction_temp;


          }
        }

        //counter_cell += 1;
      }
    }
  }

  }// end of pragma


}

void LevelSet3D::ComputeGradients()
{
  // compute gradients using HJWENO
  // http://web.stanford.edu/class/cs237c/Lecture16.pdf

  // resize gradient vector
  grid_gradient.resize(num_grid_pts, 1.0);

  #pragma omp parallel
  {

  #pragma omp for nowait
  for(int i1 = 0; i1 < num_grid_pts; i1 ++)
  {


    if(phi_status[i1] == 1)
    {
      // this point is in the narrow band

      int sign = grid_vel[i1] < 0 ? -1 : 1;

      int x = Index_to_grid_pt_zyx( i1 )[0];
      int y = Index_to_grid_pt_zyx( i1 )[1];
      int z = Index_to_grid_pt_zyx( i1 )[2];

      //x, y, z are the values of this grid point

      //Steps: calculate grad_x_plus, grad_x_minus, grad_y_plus, grad_y_minus, grad_z_plus, grad_z_minus

      // Gradient is simply norm of gradients in all three directions



      double v1, v2, v3, v4, v5;

      double grad_x_plus, grad_y_plus, grad_z_plus , grad_x_minus, grad_y_minus, grad_z_minus;

      // -------------------------- Calulate grad_x_minus here -------------------------------------------

      if (x == 0)
      {
        // point is on the left edge
        v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;

        v3 = v4;
        v2 = v4;
        v1 = v4;
      }

      if (x == 1)
      {
        // point is on one point from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;

        v2 = v3;
        v1 = v3;
      }

      if (x == 2)
      {
        // point is two points from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;

        v1 = v2;
      }

      if (x >= 3 && x <= num_elem_x - 2)
      {
        // point is in bulk

        v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 3, y , z)] ;
      }

      if (x == num_elem_x - 1)
      {
        // point is one point from right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 3, y , z)] ;

        v5 = v4;
      }

      if (x == num_elem_x )
      {
        // point is on the right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 3, y , z)] ;

        v5 = v3;
        v4 = v3;
      }

      grad_x_minus = sign*gradHJWENO(v1, v2, v3, v4, v5);


      // -------------------------- Calulate grad_x_plus here -------------------------------------------

      if (x == 0)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;

        v4 = v3;
        v5 = v3;

      }

      if (x == 1)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;

        v5 = v4;
      }

      if (x >= 2 && x <= num_elem_x - 3)
      {
        // point is in bulk

        v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;
      }

      if (x == num_elem_x - 2)
      {
        // point is two points from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;

        v1 = v2;
      }

      if (x == num_elem_x - 1)
      {
        // point is one point from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;

        v1 = v3;
        v2 = v3;
      }


      if (x == num_elem_x )
      {
        // point is on the right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x - 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x - 2, y , z)] ;

        v1 = v4;
        v2 = v4;
        v3 = v4;
      }

      grad_x_plus = sign*gradHJWENO(v1, v2, v3, v4, v5);



      // -------------------------- Calulate grad_y_minus here -------------------------------------------

      if (y == 0)
      {
        // point is on the left edge
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y + 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 1 , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 0, z)] ;

        v3 = v4;
        v2 = v4;
        v1 = v4;
      }

      if (y == 1)
      {
        // point is on one point from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y + 2 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 1 , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y + 1 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y - 1 , z)] ;

        v2 = v3;
        v1 = v3;
      }

      if (y == 2)
      {
        // point is two points from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y + 2 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 1 , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x ,y + 1,  z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y - 1 , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y - 1 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y - 2 , z)] ;

        v1 = v2;
      }

      if (y >= 3 && y <= num_elem_y - 2)
      {
        // point is in bulk

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y + 1, z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y + 1 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y - 2 , z)] -  phi_temp[Grid_pt_to_index_zyx(x , y- 3 , z)] ;
      }

      if (y == num_elem_y - 1)
      {
        // point is one point from right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y + 0, z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y - 3, z)] ;

        v5 = v4;
      }

      if (y == num_elem_y )
      {
        // point is on the right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y + 0 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x , y - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x , y- 2 , z)] -  phi_temp[Grid_pt_to_index_zyx(x , y - 3 , z)] ;

        v5 = v3;
        v4 = v3;
      }

      grad_y_minus = sign*gradHJWENO(v1, v2, v3, v4, v5);


      // -------------------------- Calulate grad_y_plus here -------------------------------------------

      if (y == 0)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y  + 3, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y+ 1 , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y+ 0 , z)] ;

        v4 = v3;
        v5 = v3;

      }

      if (y == 1)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y + 3 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y + 2 , z)] -  phi_temp[Grid_pt_to_index_zyx(x , y + 1, z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x , y+ 0 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;

        v5 = v4;
      }

      if (y >= 2 && y <= num_elem_y - 3)
      {
        // point is in bulk

        v1 = phi_temp[Grid_pt_to_index_zyx(x, y + 3 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x , y+ 1 , z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x , y - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;
      }

      if (y == num_elem_y - 2)
      {
        // point is two points from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y  + 2, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y + 0, z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;

        v1 = v2;
      }

      if (y == num_elem_y - 1)
      {
        // point is one point from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y  + 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y  + 0, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y - 1 , z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y  - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x, y  - 2, z)] ;

        v1 = v3;
        v2 = v3;
      }


      if (y == num_elem_y )
      {
        // point is on the right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x , y + 0, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y - 1, z)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x , y - 1, z)] -  phi_temp[Grid_pt_to_index_zyx(x , y - 2, z)] ;

        v1 = v4;
        v2 = v4;
        v3 = v4;
      }

      grad_y_plus = sign*gradHJWENO(v1, v2, v3, v4, v5);


      // -------------------------- Calulate grad_z_minus here -------------------------------------------

      if (z == 0)
      {
        // point is on the left edge
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y, z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y, z + 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z + 0)] ;

        v3 = v4;
        v2 = v4;
        v1 = v4;
      }

      if (z == 1)
      {
        // point is on one point from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z+ 1 )] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z+ 0)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] ;

        v2 = v3;
        v1 = v3;
      }

      if (z == 2)
      {
        // point is two points from left edge

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z+ 1)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x ,y ,  z+ 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z+ 0)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z- 2)] ;

        v1 = v2;
      }

      if (z >= 3 && z <= num_elem_z - 2)
      {
        // point is in bulk

        v5 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y, z + 1)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z + 0 )] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z  - 1)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y, z  - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 2)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y , z - 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z- 3)] ;
      }

      if (z == num_elem_z - 1)
      {
        // point is one point from right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 1)] -  phi_temp[Grid_pt_to_index_zyx(x , y, z + 0)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z  - 1)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z  - 2)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y, z  - 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y, z - 3)] ;

        v5 = v4;
      }

      if (z == num_elem_z )
      {
        // point is on the right edge

        //v5 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v4 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z  - 1)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x , y, z - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y, z  - 2)] ;
        v1 = phi_temp[Grid_pt_to_index_zyx(x , y , z- 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z - 3)] ;

        v5 = v3;
        v4 = v3;
      }

      grad_z_minus = sign*gradHJWENO(v1, v2, v3, v4, v5);


      // -------------------------- Calulate grad_z_plus here -------------------------------------------

      if (z == 0)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 3)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z+ 1)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 1)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z+ 0)] ;

        v4 = v3;
        v5 = v3;

      }

      if (z == 1)
      {
        // point is on the left edge
        v1 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 3)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z+ 1)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y, z  + 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x , y , z+ 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] ;

        v5 = v4;
      }

      if (z >= 2 && z <= num_elem_z - 3)
      {
        // point is in bulk

        v1 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 3)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x , y , z+ 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x , y , z- 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 2)] ;
      }

      if (z == num_elem_z - 2)
      {
        // point is two points from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        v2 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 2)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z+ 0)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z- 2)] ;

        v1 = v2;
      }

      if (z == num_elem_z - 1)
      {
        // point is one point from right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        v3 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z + 0)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x, y , z + 0)] -  phi_temp[Grid_pt_to_index_zyx(x, y  , z- 1)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x, y , z - 1)] -  phi_temp[Grid_pt_to_index_zyx(x, y , z - 2)] ;

        v1 = v3;
        v2 = v3;
      }


      if (z == num_elem_z )
      {
        // point is on the right edge

        //v1 = phi_temp[Grid_pt_to_index_zyx(x + 3, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] ;
        //v2 = phi_temp[Grid_pt_to_index_zyx(x + 2, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] ;
        //v3 = phi_temp[Grid_pt_to_index_zyx(x + 1, y , z)] -  phi_temp[Grid_pt_to_index_zyx(x + 0, y , z)] ;
        v4 = phi_temp[Grid_pt_to_index_zyx(x , y , z+ 0)] -  phi_temp[Grid_pt_to_index_zyx(x , y , z - 1)] ;
        v5 = phi_temp[Grid_pt_to_index_zyx(x , y , z- 1)] -  phi_temp[Grid_pt_to_index_zyx(x , y, z - 2)] ;

        v1 = v4;
        v2 = v4;
        v3 = v4;
      }

      grad_z_plus = sign*gradHJWENO(v1, v2, v3, v4, v5);

      double grad = 0.0;

      if (grad_x_minus < 0)   grad += grad_x_minus * grad_x_minus;
      if (grad_y_minus < 0)   grad += grad_y_minus * grad_y_minus;
      if (grad_z_minus < 0)   grad += grad_z_minus * grad_z_minus;

      if (grad_x_plus > 0)   grad += grad_x_plus * grad_x_plus;
      if (grad_y_plus > 0)   grad += grad_y_plus * grad_y_plus;
      if (grad_z_plus > 0)   grad += grad_z_plus * grad_z_plus;

      grad = sqrt(grad);

      if(grad < 0.001) grad = 1.0;// this is done so that corner points have appropriate gradient value

      grid_gradient[i1] = grad;

      //cout << "x, y ,z = " << x << " " << y << " " << z <<  " ; grad = " << grad << endl;



    } // end of if loop checking narrow band status

  } // end of grid point loop

}// end of pragma loop



}

double LevelSet3D::gradHJWENO(double v1, double v2, double v3, double v4, double v5)
{
    // Approximate the gradient using the 5th order Hamilton-Jacobi WENO approximation.
    // Taken from pages 34-35 of "Level Set Methods and Dynamic Implicit Surfaces".
    // See: http://web.stanford.edu/class/cs237c/Lecture16.pdf


    double oneQuarter        = 1.0  / 4.0;
    double thirteenTwelths   = 13.0 / 12.0;
    double eps               = 1e-6;

    // Estimate the smoothness of each stencil.

    double s1 = thirteenTwelths * (v1 - 2*v2 + v3)*(v1 - 2*v2 + v3)
              + oneQuarter * (v1 - 4*v2 + 3*v3)*(v1 - 4*v2 + 3*v3);

    double s2 = thirteenTwelths * (v2 - 2*v3 + v4)*(v2 - 2*v3 + v4)
              + oneQuarter * (v2 - v4)*(v2 - v4);

    double s3 = thirteenTwelths * (v3 - 2*v4 + v5)*(v3 - 2*v4 + v5)
              + oneQuarter * (3*v3 - 4*v4 + v5)*(3*v3 - 4*v4 + v5);

    // Compute the alpha values for each stencil.

    double alpha1 = 0.1 / ((s1 + eps)*(s1 + eps));
    double alpha2 = 0.6 / ((s2 + eps)*(s2 + eps));
    double alpha3 = 0.3 / ((s3 + eps)*(s3 + eps));

    // Calculate the normalised weights.

    double totalWeight = alpha1 + alpha2 + alpha3;

    double w1 = alpha1 / totalWeight;
    double w2 = alpha2 / totalWeight;
    double w3 = alpha3 / totalWeight;

    // Sum the three stencil components.
    double grad = w1 * (2*v1 - 7*v2 + 11*v3)
                + w2 * (5*v3 - v2 + 2*v4)
                + w3 * (2*v3 + 5*v4 - v5);

    grad *= (1.0 / 6.0);

    return grad;
}

void LevelSet3D::Advect()
{

  // compute gradinets
  ComputeGradients();

  #pragma omp parallel for
  for(int i = 0; i < num_grid_pts; i++)
  {
    phi[i].val = 1.0*phi_temp[i] + 1.0*grid_vel[i] * grid_gradient[i] ;
  }

  // Make sure the phi is 0 on the domain boundary if it's positive

  #pragma omp parallel for
  for (uint i = 0; i < num_elem_x + 1; i++)
  {
    for (uint j = 0; j < num_elem_y + 1; j++)
    {
      for (uint k = 0; k < num_elem_z + 1; k++)
      {
        bool is_on_domain = (i == 0 || i == num_elem_x) ||  (j == 0 || j == num_elem_y) ||   (k == 0 || k == num_elem_z);
        if(is_on_domain && phi[Grid_pt_to_index_zyx(i,j,k)].val > 0)
        {
          phi[Grid_pt_to_index_zyx(i,j,k)].val = 0.0;
        }
      }
    }
  }
}

void LevelSet3D::Advect_LBeam()
{
  for(int i = 0; i < num_grid_pts; i++)
  {
    phi[i].val = 1.0*phi_temp[i] + 1.0*grid_vel[i] ;
  }

  // Make sure the phi is 0 on the domain boundary if it's positive

  for (uint i = 0; i < num_elem_x + 1; i++)
  {
    for (uint j = 0; j < num_elem_y + 1; j++)
    {
      for (uint k = 0; k < num_elem_z + 1; k++)
      {
        bool is_on_domain = (i == 0 || i == num_elem_x) ||  (j == 0 || j == num_elem_y) ||   (k == 0 || k == num_elem_z) || (i == nx1 && k >= nz1) || (i >= nx1 && k == nz1);
        if(is_on_domain && phi[Grid_pt_to_index_zyx(i,j,k)].val > 0)
        {
          phi[Grid_pt_to_index_zyx(i,j,k)].val = 0.0;
        }
      }
    }
  }
}

int LevelSet3D::Grid_pt_to_index_zyx(int x, int y, int z)
{
  // returns the index of a grid point in the zyx system (mp4vector convention)

  return z + y*(num_elem_z + 1) + x*(num_elem_z + 1)*(num_elem_y + 1);
}

std::vector<int> LevelSet3D::Index_to_grid_pt_zyx(int index)
{
    // returns the grid point of an index in the zyx system (mp4vector convention)

    std::vector<int> ret(3,0);

    ret[0] = int(index/((num_elem_z+1)*(num_elem_y + 1)));

    ret[1] = int((index - ret[0]*(num_elem_z+1)*(num_elem_y + 1) ) / (num_elem_z+1));

    ret[2] = index - ret[0]*(num_elem_z+1)*(num_elem_y + 1) - ret[1]*(num_elem_z+1);

    return ret;
}

void LevelSet3D::SetBoxDimensions(uint box_x, uint box_y, uint box_z)
{
  num_elem_x = box_x;
  num_elem_y = box_y;
  num_elem_z = box_z;

  num_grid_pts = (box_x+1)*(box_y+1)*(box_z+1);

  num_cells =(box_x)*(box_y)*(box_z);

}

void LevelSet3D::MakeBox()
{


  /*
   phi is the signed distance function corresponding to the box.
    phi is an array of length (num_elem_x+1)*(num_elem_y+1)*(num_elem_z+1)
    of mp4Vector points containing coordinates and signed distance values.
    The phi array is arranged so that the points first start on the z axis,
    and then go to they y, and then x axis. Phi is computed using the
    minimum distance from the boundary.
   */

  uint count_phi = 0;


  // Loop through all level set grid points
  for (uint i = 0; i < num_elem_x + 1; i++)
  {
    for (uint j = 0; j < num_elem_y + 1; j++)
    {
      for (uint k = 0; k < num_elem_z + 1; k++)
      {

        // Assign the current coordinate
        phi[count_phi].x = i;
        phi[count_phi].y = j;
        phi[count_phi].z = k;

        // Calculate the minimum distance between the current grid point,
        // at coordinate (i, j, k), and all of the boundary points

        // calculating the minimum distance this way is computationally expensive
        // use a simpler method for a box
        // just measure the minimum distance from the domain

        phi[count_phi].val =   std::min({i, num_elem_x-i, j, num_elem_y-j, k, num_elem_z - k});

        // create holes here

        int hole_count = holes.size();

        for(int hole_iter = 0; hole_iter < hole_count; hole_iter++)
        {

          std::vector<double> hole =  holes[hole_iter];

          float dist = (hole[0] - i)*(hole[0] - i) + (hole[1] - j)*(hole[1] - j) + (hole[2] - k)*(hole[2] - k) - hole[3]*hole[3];

          if(dist >= 0) dist = std::sqrt(dist);
          if(dist < 0) dist = -std::sqrt(-dist);

          phi[count_phi].val = std::min(phi[count_phi].val, dist );

        }

        count_phi++;

      }
    }
  }



  std::cout << "\nNumber of phi values: " << num_grid_pts << std::endl;


}

void LevelSet3D::MakeLBeam()
{

  //uint nx1 = 16; // length of the vertical branch of the L
  //uint nz1 = 16; // length of the horizontal branch of the L

  /*
   * phi is the signed distance function corresponding to the box.
   * phi is an array of length (num_elem_x+1)*(num_elem_y+1)*(num_elem_z+1)
   * of mp4Vector points containing coordinates and signed distance values.
   * The phi array is arranged so that the points first start on the z axis,
   * and then go to they y, and then x axis. Phi is computed using the
   * minimum distance from the boundary points.
   */

  uint count_phi = 0;


  // Loop through all level set grid points
  for (uint i = 0; i < num_elem_x + 1; i++)
  {
    for (uint j = 0; j < num_elem_y + 1; j++)
    {
      for (uint k = 0; k < num_elem_z + 1; k++)
      {

        // Assign the current coordinate
        phi[count_phi].x = i;
        phi[count_phi].y = j;
        phi[count_phi].z = k;

        // Calculate the minimum distance between the current grid point,
        // at coordinate (i, j, k), and all of the boundary points

        // calculating the minimum distance this way is computationally expensive
        // use a simpler method for the L
        // just measure the minimum distance from the domain
        // divide region into 4 quadrants and measure minimum distance

        if( i <= nx1 && k <= nz1)
        {
          // this region is the bottom left hand side

          phi[count_phi].val =   std::min({i, j, num_elem_y-j, k});

          float corner_dist = std::sqrt( (nx1 - i)*(nx1 - i) + (nz1 - k)*( nz1 - k)  );

          phi[count_phi].val =   std::min( phi[count_phi].val , corner_dist);

        }
        else if (i > nx1 && k <= nz1)
        {
          // this region is the bottom right hand side

          phi[count_phi].val =   std::min({num_elem_x-i, j, num_elem_y-j, k, nz1 - k});

        }
        else if (i <= nx1 && k > nz1)
        {
          // this region is the top left side

          phi[count_phi].val =   std::min({i, nx1-i, j, num_elem_y-j, nz1 - k});
        }
        else
        {
          // this region is the top right side where the hole is located

          phi[count_phi].val =   std::min({i - nx1, k - nz1});

          phi[count_phi].val = - phi[count_phi].val;


        }

        //phi[count_phi].val =   std::min({i, num_elem_x-i, j, num_elem_y-j, k, num_elem_z - k});

        // create holes here


        int hole_count = holes.size();

        for(int hole_iter = 0; hole_iter < hole_count; hole_iter++)
        {

          std::vector<double> hole =  holes[hole_iter];

          float dist = (hole[0] - i)*(hole[0] - i) + (hole[1] - j)*(hole[1] - j) + (hole[2] - k)*(hole[2] - k) - hole[3]*hole[3];

          if(dist >= 0) dist = std::sqrt(dist);
          if(dist < 0) dist = -std::sqrt(-dist);

          phi[count_phi].val = std::min(phi[count_phi].val, dist );

        }

        count_phi++;

      }
    }
  }



  std::cout << "\nNumber of phi values: " << num_grid_pts << std::endl;


}

void LevelSet3D::MarchingCubesWrapper()
{
  /*
   * MarchingCubesCross is used to calculate the zero surface, which
   * implements the marching cubes method. It requires: the number of
   * elements in the x, y, and z directions; the value at which the
   * surface occurs, so zero in our case; an array of mp4Vector
   * points which contain the grid coordinates and associated signed
   * distance values; and a reference to a variable that will be
   * updated after the method finishes. The method returns a pointer
   * to an array of triangle elements and updates a variable with the
   * number of triangles in the array.
   */
   //cout << "Echo 1a " << endl;
   triangle_array = MarchingCubesCross(num_elem_x, num_elem_y, num_elem_z, iso_value, phi, num_triangles);
   //cout << "Echo 1b " << endl;

    /*
      It's computationally expensive to compute the vertices without repititions.
      Instead, compute midpoints and assume they are the boundary points

      Also, compute areas of the triangles here

    */


   // Resize boundary_pts
   boundary_pts.resize(0);
   boundary_areas.resize(0);
   boundary_pts_one_vector.resize(0);
   for(int i = 0; i < num_triangles; i++)
   {
     std::vector<double> point_temp(3,0);

     point_temp[0] += triangle_array[i].p[0].x;
     point_temp[0] += triangle_array[i].p[1].x;
     point_temp[0] += triangle_array[i].p[2].x;
     point_temp[0] = point_temp[0]/3.0;

     point_temp[1] += triangle_array[i].p[0].y;
     point_temp[1] += triangle_array[i].p[1].y;
     point_temp[1] += triangle_array[i].p[2].y;
     point_temp[1] = point_temp[1]/3.0;

     point_temp[2] += triangle_array[i].p[0].z;
     point_temp[2] += triangle_array[i].p[1].z;
     point_temp[2] += triangle_array[i].p[2].z;
     point_temp[2] = point_temp[2]/3.0;

     //cout << point_temp[0] << " "<< point_temp[1] <<" " << point_temp[2] << endl;

     boundary_pts.push_back(point_temp);

     boundary_pts_one_vector.push_back(point_temp[0]);
     boundary_pts_one_vector.push_back(point_temp[1]);
     boundary_pts_one_vector.push_back(point_temp[2]);

     // Calculate the area of the triangle using the good old Heron's formula

     double a = std::pow(triangle_array[i].p[0].x - triangle_array[i].p[1].x,2);
     a+=  std::pow(triangle_array[i].p[0].y - triangle_array[i].p[1].y,2);
     a+=  std::pow(triangle_array[i].p[0].z - triangle_array[i].p[1].z,2);
     a = std::sqrt(a);

     double b = std::pow(triangle_array[i].p[0].x - triangle_array[i].p[2].x,2);
     b+=  std::pow(triangle_array[i].p[0].y - triangle_array[i].p[2].y,2);
     b+=  std::pow(triangle_array[i].p[0].z - triangle_array[i].p[2].z,2);
     b = std::sqrt(b);

     double c = std::pow(triangle_array[i].p[2].x - triangle_array[i].p[1].x,2);
     c+=  std::pow(triangle_array[i].p[2].y - triangle_array[i].p[1].y,2);
     c+=  std::pow(triangle_array[i].p[2].z - triangle_array[i].p[1].z,2);
     c = std::sqrt(c);

     double s = 0.5*(a+b+c);

     double area = s*(s-a)*(s-b)*(s-c);

     area = std::sqrt(std::max(area,0.0));

     //cout << area << endl;

     boundary_areas.push_back(area);

   }

   //cout << "Boundary points size = " << boundary_pts.size() << endl;
   num_boundary_pts = boundary_pts.size();
}

void LevelSet3D::ExtrapolateVelocities()
{

  /*
    Extrapolate velocities to grid points

    Requires a set of optimum velocities along the boundary points

    Let's not be too fancy here. Just use the good old inverse square method

  */


  std::vector<float> weight(num_grid_pts,0.0);
  std::vector<float> weightedvel(num_grid_pts,0.0);



  grid_vel.resize(num_grid_pts , 0.0); // Reinitialize to zeros


  int vel_band_width = 2+0*narrow_band_width;


  #pragma omp parallel for
  for(int i1 = 0; i1 < num_boundary_pts; i1++)
  {

    for(int i = 1-vel_band_width; i <= 1+vel_band_width; i++)
    {
      for(int j = 1-vel_band_width; j <= 1+vel_band_width; j++)
      {
        for(int k = 1-vel_band_width; k <= 1+vel_band_width; k++)
        {

          std::vector<double> current_pt(3,0);
          current_pt[0] = boundary_pts[i1][0];
          current_pt[1] = boundary_pts[i1][1];
          current_pt[2] = boundary_pts[i1][2];

          /*
            x_index, y_index, z_index are grid cordinates + 1

            [x_index,y_index, z_index] is [1,1,1] for a point[0,0,0]

            Stupid Matlab!
          */



          int x_index  = floor(current_pt[0]+0.5) + i;
          int y_index  = floor(current_pt[1]+0.5) + j;
          int z_index  = floor(current_pt[2]+0.5) + k;

          if(x_index > 0 && y_index >0 && z_index > 0 &&  x_index < num_elem_x+2 &&  y_index < num_elem_y+2 &&   z_index < num_elem_z+2)
          {
            float dist = (x_index - 1 - current_pt[0])*(x_index - 1 - current_pt[0]);
            dist += (y_index - 1 - current_pt[1])*(y_index - 1 - current_pt[1]);
            dist += (z_index - 1 - current_pt[2])*(z_index - 1 - current_pt[2]);

            dist =std::sqrt(dist);

            float weight_temp = 1.0/std::max(dist,float(0.000001));
            weight_temp = weight_temp*weight_temp;

            /*
              Check if the point has alredy been assigned a weight
              If not, then assign weight and velocity
              If yes, then update weight and velocity
            */
            int index_current_pt = Grid_pt_to_index_zyx(x_index-1,y_index-1,z_index-1);

            weightedvel[index_current_pt] += weight_temp*opt_vel[i1];
            weight[index_current_pt] += weight_temp;

          }
        }
      }

    }
  }//end of i1

  #pragma omp parallel for
  for(int i = 0; i < num_grid_pts; i++){
    if(weight[i] > 0) grid_vel[i] = weightedvel[i]/weight[i];
  }


}

void LevelSet3D::SolveEikonal(std::vector<int> indices_xyz)
{

  int index_x = indices_xyz[0];
  int index_y = indices_xyz[1];
  int index_z = indices_xyz[2];



  //cout << "Eikonal indices are "<< index_x << " " << index_y << " "<< index_z << endl;
  // phix, phiy, phiz are the adjacent values of the boundary point

  double phix, phiy, phiz;

  if(index_x == 0)
  {
   phix = phi_temp[Grid_pt_to_index_zyx(index_x+1,index_y,index_z)];
  }
  else if (index_x == num_elem_x)
  {
   phix = phi_temp[Grid_pt_to_index_zyx(index_x-1,index_y,index_z)];
  }
  else
  {
   phix = std::min( phi_temp[Grid_pt_to_index_zyx(index_x+1,index_y,index_z)] , phi_temp[Grid_pt_to_index_zyx(index_x-1,index_y,index_z)] ) ;
  }

  if(index_y == 0)
  {
   phiy = phi_temp[Grid_pt_to_index_zyx(index_x,index_y+1,index_z)];
  }
  else if (index_y == num_elem_y)
  {
   phiy = phi_temp[Grid_pt_to_index_zyx(index_x,index_y-1,index_z)];
  }
  else
  {
   phiy = std::min( phi_temp[Grid_pt_to_index_zyx(index_x,index_y-1,index_z)] , phi_temp[Grid_pt_to_index_zyx(index_x,index_y+1,index_z)] ) ;
  }

  if(index_z == 0)
  {
   phiz = phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z+1)];
  }
  else if (index_z == num_elem_z)
  {
   phiz = phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z-1)];
  }
  else
  {
   phiz = std::min( phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z-1)] , phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z+1)] ) ;
  }

  double a_quad = 3.0;
  double b_quad = -2.0*(phix + phiy + phiz);
  double c_quad = phix*phix + phiy*phiy + phiz*phiz - 1.0;

  if(b_quad*b_quad >= 4*a_quad*c_quad)
  {
    phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z)] = 0.5*(-b_quad + std::sqrt(b_quad*b_quad - 4*a_quad*c_quad))/a_quad;
  }
  else
  {
    phi_temp[Grid_pt_to_index_zyx(index_x,index_y,index_z)] = std::min({phix,phiy,phiz}) + 0.75;
  }





}

void LevelSet3D::SetupNarrowBand()
{
  /*
    Set up narrow band.
      phi_status is:
      -1 if its out of the narrow band
      0 for Far
      1 for accepted
      2 for Considered

     Assume initially all points are far

  */

  phi_status.resize(0);
  phi_status.resize(num_grid_pts,-1);



  /*
    Reset phi_temp. Assign large values based on whether or not
    a point is inside or outside.
    Save phi values from previous iteration for accepted points
    And then save points as considered

  */

  //cout << "Echo 3" << endl;

  /*
    Copy phi_temp from previous iteration
  */
  phi_temp.resize(0);
  phi_temp.resize(num_grid_pts,0.0);

  //cout << "Echo 4" << endl;

  #pragma omp parallel for
  for(int i1 = 0; i1 < num_grid_pts; i1++)
  {
    phi_temp[i1] = phi[i1].val;
  }



  #pragma omp parallel
  {

  #pragma omp for nowait
  for(int i1 = 0; i1 < num_boundary_pts; i1++)
  {

    for(int i = 1-narrow_band_width; i <= 1+narrow_band_width; i++)
    {
      for(int j = 1-narrow_band_width; j <= 1+narrow_band_width; j++)
      {
        for(int k = 1-narrow_band_width; k <= 1+narrow_band_width; k++)
        {

          std::vector<double> current_pt(3,0);
          current_pt[0] = boundary_pts[i1][0];
          current_pt[1] = boundary_pts[i1][1];
          current_pt[2] = boundary_pts[i1][2];

          int x_index  = floor(current_pt[0]+0.5) + i;
          int y_index  = floor(current_pt[1]+0.5) + j;
          int z_index  = floor(current_pt[2]+0.5) + k;

          if(x_index > 0 && y_index >0 && z_index > 0 &&  x_index < num_elem_x+2 &&  y_index < num_elem_y+2 &&   z_index < num_elem_z+2)
          {
            //grid point is inside domain

            if( std::max( { std::abs(x_index - 1 - current_pt[0]),std::abs( y_index - 1 - current_pt[1]), std::abs(z_index - 1 - current_pt[2]) }  ) < 1.0001)
            {
              // save point as accepted
              phi_status[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] = 1;

              // copy phi from previous iteration
              //phi_temp[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] = phi[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)].val;
            }

          }

        }
      }
    }

  }

} // end of pragma




  //Now collect all the considered points

  // Resize indices_considered
  indices_considered_inside.resize(0);
  indices_considered_outside.resize(0);

  #pragma omp parallel
  {

    std::vector<int> indices_considered_inside_private;
    std::vector<int> indices_considered_outside_private;

    #pragma omp for nowait //fill vec_private in parallel
    for(int i1 = 0; i1 < num_boundary_pts; i1++)
    {

      for(int i = 1-narrow_band_width; i <= 1+narrow_band_width; i++)
      {
        for(int j = 1-narrow_band_width; j <= 1+narrow_band_width; j++)
        {
          for(int k = 1-narrow_band_width; k <= 1+narrow_band_width; k++)
          {

            std::vector<double> current_pt(3,0);
            current_pt[0] = boundary_pts[i1][0];
            current_pt[1] = boundary_pts[i1][1];
            current_pt[2] = boundary_pts[i1][2];

            int x_index  = floor(current_pt[0]+0.5) + i;
            int y_index  = floor(current_pt[1]+0.5) + j;
            int z_index  = floor(current_pt[2]+0.5) + k;

            if( std::max( { std::abs(x_index - 1 - current_pt[0]),std::abs( y_index - 1 - current_pt[1]), std::abs(z_index - 1 - current_pt[2]) }  ) <= narrow_band_width)
            {
              //grid point is in the narrow band

              if(x_index > 0 && y_index >0 && z_index > 0 &&  x_index < num_elem_x+2 &&  y_index < num_elem_y+2 &&   z_index < num_elem_z+2)
              {
                //grid point is inside domain

                //Copy phi values from previous  iteration
                //phi_temp[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] = phi[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)].val;


                // If not already considered or accepted
                // save indices of these "considered" points
                // set status as considererd

                if(phi_status[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] < 1)
                {
                  if( phi[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)].val >= 0)
                  {
                    phi_status[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] = 2;
                    indices_considered_inside_private.push_back(Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1));

                    //cout << "indices inside = " << x_index-1 << " " << y_index-1 <<" " << z_index-1<<  endl;
                  }
                  else
                  {
                    phi_status[Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1)] = 2;
                    indices_considered_outside_private.push_back(Grid_pt_to_index_zyx(x_index-1 , y_index - 1, z_index - 1));
                  }

                }
              }

            }



          }
        }
      }


    }

    #pragma omp critical
    indices_considered_inside.insert(indices_considered_inside.end(), indices_considered_inside_private.begin(), indices_considered_inside_private.end());

    #pragma omp critical
    indices_considered_outside.insert(indices_considered_outside.end(), indices_considered_outside_private.begin(), indices_considered_outside_private.end());



  }








}

void LevelSet3D::FastMarchingMethod()
{

  //Compute phi_considered for the considered points by solving the
  //quadratic equation


  phi_considered.resize(indices_considered.size());


  #pragma omp parallel for
  for(int i1 = 0; i1 < indices_considered.size(); i1 ++)
  {
    //solve the eikonal equation at this grid point
    SolveEikonal(Index_to_grid_pt_zyx( indices_considered[i1] ));
    phi_considered[i1] = (phi_temp[indices_considered[i1]]);


  }

  // Step-2
  // Sort phi_considered and save the indices

  std::size_t n(0);
  std::vector<int> y(indices_considered.size());
  std::generate(std::begin(y), std::end(y), [&]{ return n++; });

  std::sort(std::begin(y), std::end(y), [&](int i1, int i2) { return phi_considered[i1] < phi_considered[i2]; } );

  std::sort(std::begin(phi_considered), std::end(phi_considered));

  std::vector<int> indices_considered_dup = indices_considered;

  #pragma omp parallel for
  for(int i =0; i<  indices_considered.size(); i++ )
  {
    indices_considered[y[i]] =  indices_considered_dup[i];
  }

  // Step-3
  // Update velocity at this point

  for(int i =0; i<  indices_considered.size(); i++){
    // Update velocity of this point
    int x_index = Index_to_grid_pt_zyx( indices_considered[i] )[0];
    int y_index = Index_to_grid_pt_zyx( indices_considered[i] )[1];
    int z_index = Index_to_grid_pt_zyx( indices_considered[i] )[2];

    SolveEikonal(Index_to_grid_pt_zyx( indices_considered[i] ));
    //phi_status[indices_considered[i]] = 1;
    UpdateVelocity(x_index, y_index, z_index);

  }

  }

void LevelSet3D::UpdateVelocity(int xin, int yin, int zin)
{

  int x_index = xin;
  int y_index = yin;
  int z_index = zin;

  double ext_weight_x, ext_weight_y, ext_weight_z;
  double ext_vel_x, ext_vel_y, ext_vel_z;

  if(x_index == 0)
  {
    // This means the point to the right of the current point is accepted
    // Use its info to update velocity

    ext_weight_x = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index + 1,y_index,z_index)];
    ext_vel_x = grid_vel[Grid_pt_to_index_zyx(x_index + 1,y_index,z_index)];

  }
  else if (x_index == num_elem_x)
  {
    ext_weight_x = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index - 1,y_index,z_index)];
    ext_vel_x = grid_vel[Grid_pt_to_index_zyx(x_index - 1,y_index,z_index)];

  }
  else
  {
    if(phi_temp[Grid_pt_to_index_zyx(x_index - 1,y_index,z_index)] < phi_temp[Grid_pt_to_index_zyx(x_index + 1,y_index,z_index)] )
    {
      // This means the point to the left of the current point is accepted
      // Use its info to update velocity
      ext_weight_x = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index -1,y_index,z_index)];
      ext_vel_x = grid_vel[Grid_pt_to_index_zyx(x_index -1,y_index,z_index)];
    }
    else
    {
      ext_weight_x = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index + 1,y_index,z_index)];
      ext_vel_x = grid_vel[Grid_pt_to_index_zyx(x_index + 1,y_index,z_index)];

    }

  }

  //similary for other dimensions

  if(y_index == 0)
  {
    // This means the point to the right of the current point is accepted
    // Use its info to update velocity

    ext_weight_y = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index + 1,z_index)];
    ext_vel_y = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index + 1,z_index)];

  }
  else if (y_index == num_elem_y)
  {
    ext_weight_y = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index - 1,z_index)];
    ext_vel_y = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index - 1,z_index)];

  }
  else
  {
    if(phi_temp[Grid_pt_to_index_zyx(x_index ,y_index - 1,z_index)] < phi_temp[Grid_pt_to_index_zyx(x_index ,y_index + 1,z_index)])
    {
      // This means the point to the left of the current point is accepted
      // Use its info to update velocity
      ext_weight_y = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index - 1,z_index)];
      ext_vel_y = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index - 1,z_index)];
    }
    else
    {
      ext_weight_y = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index  ,y_index + 1,z_index)];
      ext_vel_y = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index + 1,z_index)];

    }

  }


  if(z_index == 0)
  {
    // This means the point to the right of the current point is accepted
    // Use its info to update velocity

    ext_weight_z = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index ,z_index  + 1)];
    ext_vel_z = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index ,z_index + 1)];

  }
  else if (z_index == num_elem_z)
  {
    ext_weight_z = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index ,z_index - 1)];
    ext_vel_z = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index ,z_index - 1)];

  }
  else
  {
    if(phi_temp[Grid_pt_to_index_zyx(x_index ,y_index ,z_index - 1)] < phi_temp[Grid_pt_to_index_zyx(x_index ,y_index ,z_index + 1)])
    {
      // This means the point to the left of the current point is accepted
      // Use its info to update velocity
      ext_weight_z = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index ,y_index ,z_index - 1)];
      ext_vel_z = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index ,z_index - 1)];
    }
    else
    {
      ext_weight_z = phi_temp[Grid_pt_to_index_zyx(x_index,y_index,z_index)] - phi_temp[Grid_pt_to_index_zyx(x_index  ,y_index ,z_index + 1)];
      ext_vel_z = grid_vel[Grid_pt_to_index_zyx(x_index ,y_index ,z_index + 1)];

    }



  }

  // set up bounds for ext_weights
  ext_weight_x = std::max(1.0e-6 , ext_weight_x );
  ext_weight_y = std::max(1.0e-6 , ext_weight_y );
  ext_weight_z = std::max(1.0e-6 , ext_weight_z );


  // update velocity
  double move_limit = 1.0; // CFL limit?
  double gridvel  = ext_vel_x*ext_weight_x + ext_vel_y*ext_weight_y + ext_vel_z*ext_weight_z;
  gridvel = gridvel/(ext_weight_x + ext_weight_y + ext_weight_z);

  gridvel = std::max(gridvel , -move_limit);
  gridvel = std::min(gridvel , move_limit);

  grid_vel[Grid_pt_to_index_zyx(x_index,y_index,z_index)]  = gridvel ;


}

void LevelSet3D::WriteSTL(int box_smooth)
{

  //smooth
  std::vector<float> smooth_phi(num_grid_pts , 0.0) ;

  int box_z = num_elem_z , box_y = num_elem_y , box_x = num_elem_x;

  for(int k = 0; k <= box_z; k ++)
  {
    for(int j = 0; j <= box_y; j ++)
    {
      for(int i = 0; i <= box_x; i ++)
      {

        int count_smooth = 0;
        for(int k1 = max( 0 , k -box_smooth ); k1 <= min(int(box_z), k+box_smooth); k1 ++)
         {
           for(int j1 = max( 0 , j -box_smooth ); j1 <= min(int(box_y), j+box_smooth); j1 ++)
           {
             for(int i1 = max( 0 , i -box_smooth ); i1 <= min(int(box_x), i+box_smooth); i1 ++)
             {
               smooth_phi[Grid_pt_to_index_zyx(i,j,k)] += phi[Grid_pt_to_index_zyx(i1,j1,k1)].val ;
               count_smooth++;

             }
           }
         }

         smooth_phi[Grid_pt_to_index_zyx(i,j,k)] = smooth_phi[Grid_pt_to_index_zyx(i,j,k)] / count_smooth ;

         bool is_on_domain = ( i == 0 || i == (int(box_x)) ||  j == 0 || j == (int(box_y)) ||   k == 0 || k == (int(box_z)) );

         if(is_on_domain && smooth_phi[Grid_pt_to_index_zyx(i,j,k)]  > 0)
         {
           smooth_phi[Grid_pt_to_index_zyx(i,j,k)] = 0.0;
         }

      }
    }

  }

  // replace phi with smooth phi
  for(int k = 0; k <= box_z; k ++)
  {
    for(int j = 0; j <= box_y; j ++)
    {
      for(int i = 0; i <= box_x; i ++)
      {
        phi[Grid_pt_to_index_zyx(i,j,k)].val = smooth_phi[Grid_pt_to_index_zyx(i,j,k)] ;
      }
    }
  }

  // Redo marching cubes
  delete[] triangle_array;
  MarchingCubesWrapper();

  // Writing stl file

  std::ofstream txtfile;

  txtfile.open("mystlfile.stl");

  txtfile << "solid mysolid" << endl;

  for(int i = 0; i <  num_triangles; i++)
  {
    txtfile << "facet normal 0 0 0" << endl;
    txtfile << "  outer loop" << endl;
    txtfile << "    vertex "<<  triangle_array[i].p[0].x << " " <<  triangle_array[i].p[0].y << " " <<  triangle_array[i].p[0].z << " " << endl;
    txtfile << "    vertex "<<  triangle_array[i].p[1].x << " " <<  triangle_array[i].p[1].y << " " <<  triangle_array[i].p[1].z << " " << endl;
    txtfile << "    vertex "<<  triangle_array[i].p[2].x << " " <<  triangle_array[i].p[2].y << " " <<  triangle_array[i].p[2].z << " " << endl;
    txtfile << "  end loop" << endl;
    txtfile << "endfacet" << endl;
  }
  txtfile << "endsolid mysolid" << endl;
  txtfile.close();

}

// ----

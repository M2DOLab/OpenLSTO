


void PerformOptimization ( SensitivityData &SensData) {

  // number of elements in x, y, and z directions
  double nx = SensData.nx;
  double ny = SensData.ny;
  double nz = SensData.nz;


  // parameters for using simplex
  std::vector<double> LB = SensData.LB;
  std::vector<double> UB = SensData.UB;
  std::vector<double> LAM = SensData.LAM;

  // maximum volume of structure, in percentage
  double MaxVol = SensData.MaxVol;

  // Copy volumefraction_vector
  std::vector<double> volumefraction_vector(int(nx*ny*nz), 1.0E-5);

  volumefraction_vector = SensData.volumeFractions;
  double Vol = 0;
  for (int i = 0 ; i < volumefraction_vector.size() ; ++i) {
     Vol+= volumefraction_vector[i];
  }
  SensData.Vol = Vol;

  // Read iteration here
  int ITER = 1;
  ITER = SensData.iter;

  // Read compliance here
  double compliance = 1.0;
  compliance = SensData.compliance;


  // Read boundary points here
  std::vector< std::vector<double> > bpoints;
  double bpoint_tempx = 0.0, bpoint_tempy = 0.0, bpoint_tempz = 0.0;
  vector<double> bPoint (3, 0) ;
  vector<double> bPointArr;

  bPointArr.resize(0);
  bPointArr = SensData.bPoints;


  // Read Boundary Sensitivites
  // bsens - > vector of compliance sensitivites
  // vsens - > vector of volume sensitivites, typically = -1
  int bpointsize = bPointArr.size()/3;
  SensData.bpointsize = bpointsize;
  std::vector<double> bsens;
  bsens = SensData.bsens;

  std::vector<double> vsens;
  vsens = SensData.vsens;


  // Read Area Vector
  std::vector<double> areavector;
  areavector.resize(0);
  areavector = SensData.pointAreas;


  std::vector<double> sf; // speed of the objective
  std::vector<double> cf; // weight of the objective
  std::vector<double> sg; // speed of the constraint
  std::vector<double> cg; // weight of the constraint

  // scale sensitivities wrt the absolute values

  double abssens = std::abs(bsens[0]);

  std::vector<double> curpt = {0,0,0};

  for (int i = 0 ; i < bpointsize ; ++i)
  {
    curpt[0] = bPointArr[3*i];
    curpt[1] = bPointArr[3*i+1];
    curpt[2] = bPointArr[3*i+2];


    sf.push_back(-bsens[i]);
    sg.push_back(vsens[i]);

    if (std::abs(bsens[i]) > abssens)
    {
       abssens = abs(bsens[i]);
    }
  }

  for (int i = 0 ; i < bpointsize ; ++i)
  {
    //sf[i] = sf[i]/abssens;
    cf.push_back(sf[i]*areavector[i]);
    cg.push_back(sg[i]*areavector[i]);

  }

  cout<< "Iter: "<< ITER<< " Obj: " <<compliance<< " Constr: " << Vol/nx/ny/nz << endl;


  if(SensData.algo == 0)
  {
    double move_limit = SensData.move_limit;

    double lambda_g = move_limit;

    double percent_vol = 0.5; // delta_area will be 50 percent of maximum area

    double target_vol = Vol;
    for(int i = 0; i<cg.size(); i++){
      target_vol += cg[i]*percent_vol*(lambda_g);
    }

    target_vol = std::max(MaxVol*nx*ny*nz/100.0, target_vol);

    // here we solve this equation for lamda_f:
    // target_area = boundary_area + sum( cg_i(lambda_f*sf_i + lambda_g*sg_i ) )


    // compute domain distance (distance from the FE domain)

    std::vector<double> domain_distance_vector(bpointsize,0.0);
    for (int i = 0; i < bpointsize; i++)
    {
      curpt[0] = bPointArr[3*i];
      curpt[1] = bPointArr[3*i+1];
      curpt[2] = bPointArr[3*i+2];

      double domdist; // distance from fem domain

      domdist = std::min({ std::abs( curpt[0] - 0.0),std::abs( curpt[0] - 1.0*nx),std::abs( curpt[1] - 0.0),std::abs( curpt[1] - 1.0*ny),std::abs( curpt[2] - 0.0),std::abs( curpt[2] - 1.0*nz)  });

      if(  ( curpt[0] - nx >=0 || -(curpt[0] -0) >=0  || curpt[1] - ny >=0 || -(curpt[1] -0) >=0  || curpt[2] - nz >=0 || -(curpt[2] -0) >=0) )
      {
        domdist = -1.0*domdist;
      }
      domain_distance_vector[i] = std::min(domdist , SensData.move_limit );

    }

    double lambda_0 = 0.0;
    double delta_lambda = 0.01;
    double lambda_cur;
    double new_vol;

    double default_vol = Vol;
    for(int i = 0; i<cg.size(); i++){
      default_vol += cg[i]*lambda_g ;
    }

    for (int iter_NR = 0; iter_NR < 50; iter_NR ++)
    {
      lambda_cur = lambda_0 + 0*delta_lambda;
      // compute new area
      double new_vol0 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol0 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      lambda_cur = lambda_0 + delta_lambda;
      // compute new area
      double new_vol2 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol2 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      lambda_cur = lambda_0 - delta_lambda;
      // compute new area
      double new_vol1 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol1 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      double slope = (new_vol2 - new_vol1)/ 2 / delta_lambda;

      lambda_0 -= (new_vol0 - target_vol)/slope;



      if(std::abs(new_vol0 - target_vol) < 1.0E-5) break;
    }



    double lambda_f = lambda_0;

    // Assign optimum velocities here

    std::vector<double> OptVel(bpointsize);


    SensData.opt_vel_nlopt.resize(bpointsize);

    double abs_vel = 0.0;

    for (int i = 0; i < bpointsize; i++)
    {
      SensData.opt_vel_nlopt[i] = std::min(lambda_f*sf[i] + lambda_g*sg[i],domain_distance_vector[i]);

      // Cap the optimum velocities between +- move limit
      SensData.opt_vel_nlopt[i] = std::min(move_limit , SensData.opt_vel_nlopt[i] ) ;
      SensData.opt_vel_nlopt[i] = std::max(-move_limit , SensData.opt_vel_nlopt[i] ) ;

      abs_vel = std::max( abs_vel, std::abs(SensData.opt_vel_nlopt[i]) );
    }


    // Implement cfl conditon here

    if (abs_vel > move_limit)
    {
      for (int i = 0; i < bpointsize; i++)
      {
        SensData.opt_vel_nlopt[i] = move_limit * SensData.opt_vel_nlopt[i] / abs_vel;
      }
    }

  }

  else
  {

    // Implement your favorite optimzer here

  }



}

void PerformOptimization_LBeam ( SensitivityData &SensData) {

  // number of elements in x, y, and z directions
  double nx = SensData.nx;
  double ny = SensData.ny;
  double nz = SensData.nz;


  // parameters for using simplex
  std::vector<double> LB = SensData.LB;
  std::vector<double> UB = SensData.UB;
  std::vector<double> LAM = SensData.LAM;

  // maximum volume of structure, in percentage
  double MaxVol = SensData.MaxVol;

  // Copy volumefraction_vector
  std::vector<double> volumefraction_vector(int(nx*ny*nz), 1.0E-5);

  volumefraction_vector = SensData.volumeFractions;
  double Vol = 0;
  for (int i = 0 ; i < volumefraction_vector.size() ; ++i) {
     Vol+= volumefraction_vector[i];
  }
  SensData.Vol = Vol;

  // Read iteration here
  int ITER = 1;
  ITER = SensData.iter;

  // Read compliance here
  double compliance = 1.0;
  compliance = SensData.compliance;


  // Read boundary points here
  std::vector< std::vector<double> > bpoints;
  double bpoint_tempx = 0.0, bpoint_tempy = 0.0, bpoint_tempz = 0.0;
  vector<double> bPoint (3, 0) ;
  vector<double> bPointArr;

  bPointArr.resize(0);
  bPointArr = SensData.bPoints;

  // set up nx1 and nz1
  int nx1 = SensData.nx1;
  int nz1 = SensData.nz1;


  // Read Boundary Sensitivites
  // bsens - > vector of compliance sensitivites
  // vsens - > vector of volume sensitivites, typically = -1
  int bpointsize = bPointArr.size()/3;
  SensData.bpointsize = bpointsize;
  std::vector<double> bsens;
  bsens = SensData.bsens;

  std::vector<double> vsens;
  vsens = SensData.vsens;


  // Read Area Vector
  std::vector<double> areavector;
  areavector.resize(0);
  areavector = SensData.pointAreas;


  std::vector<double> sf; // speed of the objective
  std::vector<double> cf; // weight of the objective
  std::vector<double> sg; // speed of the constraint
  std::vector<double> cg; // weight of the constraint

  // scale sensitivities wrt the absolute values

  double abssens = std::abs(bsens[0]);

  std::vector<double> curpt = {0,0,0};

  for (int i = 0 ; i < bpointsize ; ++i)
  {
    curpt[0] = bPointArr[3*i];
    curpt[1] = bPointArr[3*i+1];
    curpt[2] = bPointArr[3*i+2];


    sf.push_back(-bsens[i]);
    sg.push_back(vsens[i]);

    if (std::abs(bsens[i]) > abssens)
    {
       abssens = std::abs(bsens[i]);
    }
  }

  for (int i = 0 ; i < bpointsize ; ++i)
  {
    sf[i] = sf[i]/abssens;
    cf.push_back(sf[i]*areavector[i]);
    cg.push_back(sg[i]*areavector[i]);

  }


  if(SensData.algo == 0)
  {
    double move_limit = SensData.move_limit;

    double lambda_g = move_limit;

    double percent_vol = 0.5; // delta_area will be 50 percent of maximum area

    double target_vol = Vol;
    for(int i = 0; i<cg.size(); i++){
      target_vol += cg[i]*percent_vol*(lambda_g);
    }

    target_vol = std::max(MaxVol*nx*ny*nz/100.0, target_vol);



    // here we solve this equation for lamda_f:
    // target_area = boundary_area + sum( cg_i(lambda_f*sf_i + lambda_g*sg_i ) )


    // compute domain distance (distance from the FE domain)

    std::vector<double> domain_distance_vector(bpointsize,0.0);
    for (int i = 0; i < bpointsize; i++)
    {
      curpt[0] = bPointArr[3*i];
      curpt[1] = bPointArr[3*i+1];
      curpt[2] = bPointArr[3*i+2];

      double domdist; // distance from fem domain

      domdist = std::min({ std::abs( curpt[0] - 0.0),std::abs( curpt[0] - nx),std::abs( curpt[1] - 0.0),std::abs( curpt[1] - ny),std::abs( curpt[2] - 0.0),std::abs( curpt[2] - nz)  });

      {
        // this region is the top right quadrant
        double domdista = std::min({( curpt[0] - 1.0*nx1),( curpt[2] - 1.0*nz1)  });

        domdista = -1.0*domdista;

        domdist = std::min(domdist , domdista);

      }

      //domdist = std::min({ std::abs( curpt[0] - 0.0),std::abs( curpt[0] - nx),std::abs( curpt[1] - 0.0),std::abs( curpt[1] - ny),std::abs( curpt[2] - 0.0),std::abs( curpt[2] - nz)  });

      domain_distance_vector[i] = domdist;

      if( std::abs(curpt[0] - nx1) + std::abs(curpt[2] - nz1) < 2.1)
      {
        //cout << " corner: " << curpt[0] << " , " << curpt[2] << "; domdist = " << domdist << endl;
      }


      //cout << " pt: " << curpt[0] << " , " <<  curpt[1] << " , "  << curpt[2] << "; domdist = " << domdist << endl;


    }

    double lambda_0 = 0.0;
    double delta_lambda = 0.01;
    double lambda_cur;
    double new_vol;

    double default_vol = Vol;
    for(int i = 0; i<cg.size(); i++){
      default_vol += cg[i]*lambda_g ;
    }

    for (int iter_NR = 0; iter_NR < 50; iter_NR ++)
    {
      lambda_cur = lambda_0 + 0*delta_lambda;
      // compute new area
      double new_vol0 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol0 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      lambda_cur = lambda_0 + delta_lambda;
      // compute new area
      double new_vol2 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol2 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      lambda_cur = lambda_0 - delta_lambda;
      // compute new area
      double new_vol1 = Vol;
      for(int i = 0; i<cg.size(); i++){
        new_vol1 -= cg[i]*std::min( domain_distance_vector[i], lambda_g*sg[i] + lambda_cur*sf[i] );
      }

      double slope = (new_vol2 - new_vol1)/ 2 / delta_lambda;

      lambda_0 -= (new_vol0 - target_vol)/slope;



      if(std::abs(new_vol0 - target_vol) < 1.0E-5) break;
    }



    double lambda_f = lambda_0;

    // Assign optimum velocities here

    std::vector<double> OptVel(bpointsize);


    SensData.opt_vel_nlopt.resize(bpointsize);

    double abs_vel = 0.0;

    for (int i = 0; i < bpointsize; i++)
    {
      SensData.opt_vel_nlopt[i] = std::min(lambda_f*sf[i] + lambda_g*sg[i],domain_distance_vector[i]);
      abs_vel = std::max( abs_vel, SensData.opt_vel_nlopt[i] );
    }

    // Implement cfl conditon here

    if (abs_vel > move_limit)
    {
      for (int i = 0; i < bpointsize; i++)
      {
        SensData.opt_vel_nlopt[i] = move_limit * SensData.opt_vel_nlopt[i] / abs_vel;
      }
    }

  }


  else
  {

    // Implement your favorite optimzer here

  }



}

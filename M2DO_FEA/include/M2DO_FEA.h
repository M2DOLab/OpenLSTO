#ifndef M2DO_FEA_MODULE_H
#define M2DO_FEA_MODULE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include <chrono>
#include <cassert>
#include <set>

#include <../../vendor/eigen3/Eigen/Dense>
#include <../../vendor/eigen3/Eigen/Sparse>

using namespace std ;
using namespace Eigen ;

#include "quadrature.h"
#include "linear_shape_function.h"
#include "node.h"
#include "heaviside_function.h"
#include "solid_element.h"
#include "solid_material.h"
#include "mesh.h"
#include "boundary_conditions.h"
#include "stationary_study.h"
#include "sensitivity.h"

#endif

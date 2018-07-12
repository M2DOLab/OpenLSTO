#ifndef M2DO_LSM_MODULE_H
#define M2DO_LSM_MODULE_H

#include <random>
#include <algorithm>
#include <functional>
#include <limits>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
//#include <nlopt.hpp>

// Should put these in the namespace!
#include "debug.h"
#include "min_unit.h"

namespace M2DO_LSM {

	#include "common.h"
	#include "mesh.h"
	#include "hole.h"
	#include "heap.h"
	#include "mersenne_twister.h"
	#include "fast_marching_method.h"
	#include "level_set.h"
	#include "boundary.h"
	#include "input_output.h"
	#include "optimise.h"
	#include "sensitivity.h"
	
}

#endif

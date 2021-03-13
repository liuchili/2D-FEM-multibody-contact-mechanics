// -*- C++ -*-
#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#include <cmath>
using std::sqrt;
using std::sin;
using std::cos;
using std::isnan;
using std::pow;
using std::cbrt;
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>		// file io
using std::stringstream;
using std::ifstream;
using std::getline;
using std::istringstream;
#include <iomanip>
#include <cassert>
#include <iostream>
#include <set>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>
#include <float.h>		// DBL_MAX, INT_MAX
#include <math.h>

#include <numeric>
using std::accumulate;
using std::cout;
using std::endl;

#include <algorithm>		// math stuff
using std::min;
using std::max;
using std::find;
using std::distance;

#include <limits>
#include <cstddef>
using std::numeric_limits;


#include <stdlib.h>     // atoi

#include <vector>  // standard vector
#include <array>
using std::vector;
using std::array;

//#include<tr1/array>
//using std::tr1::array;
using std::string;

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

#include<Eigen/SparseLU>
#include<Eigen/SparseQR>
#include <Eigen/OrderingMethods>
using Eigen::SparseLU;
using Eigen::SparseQR;
using Eigen::SparseMatrix;
using Eigen::Triplet;
#endif // DEFINITIONS_H

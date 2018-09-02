/*
 * TS: FEM GPU COMP
 * Written By: Robert Masti
 * 08/22/2018
 */
#ifndef mhdRT
#define mhdRT

#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cstring>

#include "stdio.h"
#include "mpi.h"

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// DEFINE MAXMIN FUNCTION

#define mymax(a,b) ((a>b)?a:b)
#define mymin(a,b) ((a<b)?a:b)

// DEFINE CONSTANTS

#define R 287.0
#define GAMMA 1.4
#define MU0 1.0
#define PI 3.14159265359
#define NEQ 8
#define DEBUG false

// DEFINE ELEMENT IDENTIFIERS

#define rhoid  0           // V1 and U1
#define uid 1              // V2
#define vid 2              // V3
#define wid 3              // V3
#define pid 4              // V4
#define bxid 5             // V5
#define byid 6             // V5
#define bzid 7             // V5
#define RKORDER 4         //rkorder


// Create the Class with Eigen mapp C arr
class Map2Eigen
{
  public:
  
  Map2Eigen(int ni, int nj, int nequ);

  double *Q_raw;

  Map<MatrixXd, 0, Stride<Dynamic, Dynamic> > Q[NEQ];
};


// CREATE STRUCTURE DEFINITION
struct constants
{
  int f_limiter;           // Choose which limiter you want
  int f_mesh;              // Which mesh to use
  int f_case;              // Which case to run
  int num_ghost;           // Number of ghost layers
  double cfl;              // Number of ghost layers
  double g;                // gravity
  int nmax;                // Maximum interation number
  int wint;                // Write Interval
  int pint;                // Print Interval
  int nx_c;                // number of cells x
  int ny_c;                // number of cells y
};

#endif

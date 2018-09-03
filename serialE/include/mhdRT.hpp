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
//#include "mpi.h"

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
  int num_ghost;           // Number of ghost layers
  double cfl;              // Number of ghost layers
  int nmax;                // Maximum interation number
  int wint;                // Write Interval
  int pint;                // Print Interval
};


// Function Prototypes
void inputMesh(MatrixXd& xn, MatrixXd& yn, MatrixXd& xc, MatrixXd& yc,const string mesh);

void extrapCopyCoords(MatrixXd& xc_g, MatrixXd& yc_g, const MatrixXd& xc, const MatrixXd& yc, constants C);

void computeArea(MatrixXd& Ai, MatrixXd& Aj, const MatrixXd& xn, const MatrixXd& yn);

void computeNormalVectors(MatrixXd& nix, MatrixXd& niy, MatrixXd& njx, MatrixXd& njy, const MatrixXd& xn, const MatrixXd& yn, const MatrixXd& Ai, const MatrixXd& Aj);

void computeVolume(MatrixXd& Volume, const MatrixXd& xn, const MatrixXd& yn);


void outputArray( string Address,  string FileName,  MatrixXd& out,  int n);

void outputArrayMap(string Address, string FileName, const Map2Eigen* out, int n); 

void initialize(Map2Eigen* V, const MatrixXd& xc_g, const MatrixXd& yc_g, constants C);

void setBC(Map2Eigen* V, const MatrixXd& nix, const MatrixXd& niy, const MatrixXd& njx, const MatrixXd& njy, MatrixXd& T, constants C);

void slipwallBC(Map2Eigen* V, const int Begin[], const int End[], const MatrixXd& nix, const MatrixXd& niy, const MatrixXd& njx, const MatrixXd& njy, MatrixXd& T, constants C);


void primToCons(Map2Eigen* U, const Map2Eigen* V);
void consToPrim(Map2Eigen* V, const Map2Eigen* U);
void cons2Prim(double V[],double U[],int size);

double computeTimeStep(MatrixXd& Volume, MatrixXd& Ai, MatrixXd& Aj, MatrixXd& nix,  MatrixXd& niy, MatrixXd& njx,  MatrixXd& njy, Map2Eigen* V, constants C);


double computeMaxSpeed(double V[]); 


void MUSCL(Map2Eigen* U_L, Map2Eigen* U_R, Map2Eigen* U_B, Map2Eigen* U_T, const Map2Eigen* U, constants C);


double SIGN(double a, double b);

double limiter(double& r, constants C);


void compute2dFlux(Map2Eigen* F, Map2Eigen* G,  Map2Eigen* U_L , Map2Eigen* U_R, Map2Eigen* U_B, Map2Eigen* U_T, MatrixXd& njx, MatrixXd& njy, MatrixXd& nix, MatrixXd& niy, constants C);

void computeFlux(double F[], double UA[], double UB[], double& nxhat, double& nyhat, int ForG);

void fFlux(double F[], double U[]);

void gFlux(double G[], double U[]);
#endif

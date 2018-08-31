/*
 * TS: FEM GPU COMP
 * Written By: Robert Masti
 * 08/22/2018
 */
#ifndef hw4_H_
#define hw4_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
//#include "math.h"
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstring>

//#include <Eigen/Dense>

using namespace std;
//using namespace Eigen;

// DEFINE MAXMIN FUNCTION

#define mymax(a,b) ((a>b)?a:b)
#define mymin(a,b) ((a<b)?a:b)

// DEFINE CONSTANTS

#define R 287.0
#define GAMMA 1.4
#define MU0 1.0
//#define PI 3.14159265359
#define PI 3.1415926
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

/* CREATE STRUCTURE DEFINITION
struct constants
{
  // Specify in main
  int f_case;              // 1 for curvelinear MMS, 2 for 30 degree inlet, 3 for NACA airfoil
  int f_mesh;              // 1-4 for all cases except curvelinear which has 1-6
  int f_supersonic;        // 1 for subsonic, 2 for supersonic
  int f_AOA;               // 1 for 0 deg, 2 for 8 deg. (only applicable for case_flag 3)
  int f_upwind;            // 1 for Van Leer FVS, and 2 for Roe FDS
  int f_limiter;           // 0 for limiter off, 1 for limiter on, 2 ....
  int f_eps;               // Epsilon for upwind 0 for 1st order, 1 for 2nd order
  double f_kap;
  int rk_order;            // Runge Kutta order (2nd or 4th order)
  int nmax;                // Maximum iteration number 1e6 typical
  int wint;                // Write Interval
  int pint;                // Print Interval
  int num_ghost;           // number of ghost cells
  int localdt;            // Use local time stepping? or false use global time stepping
  double tol;              // Residual tolerance 1e-10 typical
  double cfl;              // CFL number used globally
};
*/

////////////// FUNCTION PROTOTYPES ////////////////////o

int main (int argc, char * argv[]);

/*
constants loadInputFile(string FileName);

double searchInputFile(string FileName, string Var);

string buildCaseFolder(constants C);

void runCase(constants C);

//void inputMesh(vector<double>& xn, vector<double>& yn, vector<double>& zn, constants C);
*/
void meshSize(int* nx_i, int* ny_i, constants C);

string readMeshName(constants C);

void getCoord(double xn[], double yn[], double xc[], double yc[], constants C);

void computeRTheta(double rc[], double thetc[], double xc[], double yc[], constants C);

void extrapCopyCoords(double xl_g[], double xr_g[], double xb_g[], double xt_g[], double yl_g[], double yr_g[], double yb_g[], double yt_g[], double xc[], double yc[], constants C); 


void computeAreasAndNormalVectors(double njx[], double njy[], double nix[], double niy[], double Aj[],  double Ai[], double xn[], double yn[], constants C);

void computeVolume(double volume[], double xn[], double yn[], constants C);

void outputArray(string Address, string FileName, double out[], int size,  int n);

void initialize(double V[], double rc[], double thetc[], constants C);

void prim2Cons(double U[], double V[], int size);

void cons2Prim(double V[], double U[], int size);

void setBC(double Vl_g[], double Vr_g[], double Vt_g[], double Vb_g[], double njx[],double njy[], double nix[], double niy[], double V[], constants C);

double computeMaxSpeed(double V[]);

double SIGN(double a, double b);

double computeTimeStep(double volume[], double Aj[], double Ai[], double njx[], double njy[], double nix[], double niy[], double V[], constants C);

void MUSCL(double lout[], double rout[], double bout[], double tout[], double Ul_g[], double Ur_g[], double Ub_g[], double Ut_g[], double U[], constants C);

double limiter(double& r, constants C);

void getThetaExtrap(double& theta_A, double& theta_B, double U5[], int& ind, int& eq, constants& C);

void computeSource(double S[], double U[], double thetc[], constants C);


void compute2dFlux(double F[], double G[], double Ul[], double Ur[], double Ub[], double Ut[], double njx[], double njy[], double nix[], double niy[], constants C);

void computeFlux(double F[], double UA[], double UB[], double& nxhat, double& nyhat, int ForG);
#endif









/*
constants loadInputFile(string FileName);

double searchInputFile(string FileName, string Var);

string buildCaseFolder(constants C);

void runCase(constants C);

//void inputMesh(vector<double>& xn, vector<double>& yn, vector<double>& zn, constants C);
*/

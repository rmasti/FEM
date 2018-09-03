/*
 * Main file refer to mhdRT_f for functions
 * TS: FEM and GPU
 * Written By: Robert Mast
 * 08/28/2018
 */


# include "mhdRT.hpp"

int main (int argc, char * argv[])
{

  constants C;
  C.f_limiter = 1;
  C.num_ghost = 3;
  C.cfl = 0.5;
  C.nmax = 10000;
  C.wint = 10;
  C.pint = 1;
  double A = (2-1)/(1.0+2.0);
  double tend = 6.0/sqrt(A*0.1*2);

  string outputFolder =  "./output/";
  string mesh = "mesh/debugMatlab.msh";

  ///////////////// READ IN THE MESH ///////////////////
  MatrixXd xn, yn; // coord for the nodal points of cell
  MatrixXd xc, yc; // coord for the cell centers requires interpolation
  inputMesh(xn, yn, xc, yc, mesh);

  // extrap to ghost cells
  MatrixXd xc_g(xc.rows()+2*C.num_ghost, xc.cols()+2*C.num_ghost);
  MatrixXd yc_g(yc.rows()+2*C.num_ghost, yc.cols()+2*C.num_ghost);
  extrapCopyCoords(xc_g, yc_g, xc, yc, C); // get ghost cell coord

  ////////////////////// SETUP VARIABLES //////////////////////
  // SETUP AREA's and Vols
  //i or x dir areas the columns of the matrix
  MatrixXd Ai(xc.rows(), xn.cols()); // sets (row,col) (j,i)
  //j or y dir areas the rows or the matrix
  MatrixXd Aj(xn.rows(), xc.cols()); // sets (row,col) (j,i)
  MatrixXd Volume(xc.rows(), xc.cols());

  // SETUP NORMAL VECS
  //Ai normal vector has x and y components (so same size as Ai)
  MatrixXd nix(Ai.rows(), Ai.cols());
  MatrixXd niy(Ai.rows(), Ai.cols());
  //Aj normal vector has x and y componenets
  MatrixXd njx(Aj.rows(), Aj.cols());
  MatrixXd njy(Aj.rows(), Aj.cols());
 
  // time variables
  MatrixXd time(1, C.nmax);
  double dt;

  int ncols = xc.cols()+2*C.num_ghost;
  int nrows = xc.rows()+2*C.num_ghost;
  


  // Direction independent
  //
  //
  // YANG YANG YANG YANG YANG THE 2 LINES BELOW
  //
  //
  Map2Eigen *U    = new Map2Eigen(nrows, ncols, NEQ);
  cout << sizeof U->Q_raw << sizeof(double) << " " << nrows*ncols*NEQ << endl;
  Map2Eigen *U_RK = new Map2Eigen(nrows, ncols, NEQ);
  Map2Eigen *V    = new Map2Eigen(nrows, ncols, NEQ);

  Map2Eigen *S    = new Map2Eigen(xc.rows(), xc.cols(), NEQ);
  
  Map2Eigen *Res  = new Map2Eigen(xc.rows(), xc.cols(), NEQ);
  MatrixXd T(xc_g.rows(), xc_g.cols()); // temperature

  //Direction dependent
  Map2Eigen *F    = new Map2Eigen(xc.rows(), xn.cols(), NEQ);
  Map2Eigen *U_L  = new Map2Eigen(xc.rows(), xn.cols(), NEQ);
  Map2Eigen *U_R  = new Map2Eigen(xc.rows(), xn.cols(), NEQ);

  Map2Eigen *G    = new Map2Eigen(xn.rows(), xc.cols(), NEQ);
  Map2Eigen *U_T  = new Map2Eigen(xn.rows(), xc.cols(), NEQ);
  Map2Eigen *U_B  = new Map2Eigen(xn.rows(), xc.cols(), NEQ);


  ////////////////// INITIALIZE //////////////////
  computeArea(Ai, Aj, xn, yn); // take nodal coords extract interface A
  computeNormalVectors(nix, niy, njx, njy,
      xn, yn, Ai, Aj); // grab the norm vec 4 computational domain
  computeVolume(Volume, xn, yn);

  cout << "Initializing" << endl;
  initialize(V, xc_g, yc_g, C);
  T = V->Q[pid].cwiseProduct(V->Q[rhoid].cwiseInverse())/R; //computeTemperature
  setBC(V, nix, niy, njx, njy, T, C);
  primToCons(U, V);

  dt = computeTimeStep(Volume, Ai, Aj, nix, niy, njx, njy, V, C);

  outputArrayMap(outputFolder, "V", V, 0);
  outputArrayMap(outputFolder, "U", U, 0);
  outputArray(outputFolder, "xc_g", xc_g, 0);
  outputArray(outputFolder, "yc_g", yc_g, 0);
  outputArray(outputFolder, "xc", xc, 0);
  outputArray(outputFolder, "yc", yc, 0);

 // cout << U_B->Q[rhoid] << endl;
  MUSCL(U_L, U_R, U_B, U_T, U, C);



  outputArrayMap(outputFolder, "UT", U_T, 0);
  outputArrayMap(outputFolder, "UR", U_R, 0);
  compute2dFlux(F, G, U_L, U_R, U_B, U_T, njx, njy, nix, niy, C);

 
  outputArrayMap(outputFolder, "F", F, 0);
  cout << "Made it here" << endl;

   return 0;
}

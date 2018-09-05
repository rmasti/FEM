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
  C.f_limiter = 7;
  C.num_ghost = 3;
  C.cfl = 0.25;
  C.nmax = 10000;
  C.wint = 100;
  C.pint = 10;
  double A = (2-1)/(1.0+2.0);
  double tend = 6.0/sqrt(A*0.1*2);

  string outputFolder =  "./output/";
  //string outputFolder =  "/mnt/c/Users/rlm78/Desktop/matlabstuff/";
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


  MatrixXd T(nrows, ncols); // temperature
  // Direction independent
  //
  //
  // YANG YANG YANG YANG YANG THE 2 LINES BELOW
  //
  //
  Map2Eigen *U    = new Map2Eigen(nrows, ncols, NEQ);
  Map2Eigen *U_RK = new Map2Eigen(nrows, ncols, NEQ);
  Map2Eigen *V    = new Map2Eigen(nrows, ncols, NEQ);

  int njc = xc.rows();
  int nic = xc.cols();
  Map2Eigen *S    = new Map2Eigen(njc, nic, NEQ);
  Map2Eigen *Res  = new Map2Eigen(njc, nic, NEQ);
  //Direction dependent
  Map2Eigen *F    = new Map2Eigen(njc, nic+1, NEQ);
  Map2Eigen *U_L  = new Map2Eigen(njc, nic+1, NEQ);
  Map2Eigen *U_R  = new Map2Eigen(njc, nic+1, NEQ);

  Map2Eigen *G    = new Map2Eigen(njc+1, nic, NEQ);
  Map2Eigen *U_T  = new Map2Eigen(njc+1, nic, NEQ);
  Map2Eigen *U_B  = new Map2Eigen(njc+1, nic, NEQ);


  ////////////////// INITIALIZE //////////////////
  computeArea(Ai, Aj, xn, yn); // take nodal coords extract interface A // Checked
  computeNormalVectors(nix, niy, njx, njy,
      xn, yn, Ai, Aj); // grab the norm vec 4 computational domain // Checked

  computeVolume(Volume, xn, yn); // Checked

  cout << "Initializing" << endl;
  initialize(V, xc_g, yc_g, C); // Checked
  //T = V->Q[pid].cwiseProduct(V->Q[rhoid].cwiseInverse())/R; //computeTemperature
  outputArrayMap(outputFolder, "V", V, -2);
  outputArrayMap(outputFolder, "U", U, -2);
 
  primToCons(U, V);
  outputArrayMap(outputFolder, "V", V, -1);
  outputArrayMap(outputFolder, "U", U, -1);
 
  setBConU(U, nix, niy, njx, njy, C);
  //setBC(V, nix, niy, njx, njy, T, C);
  //primToCons(U, V);
  U_RK = U;

  //dt = computeTimeStep(Volume, Ai, Aj, nix, niy, njx, njy, V, C);
  dt = computeTimeStepU(Volume, Ai, Aj, nix, niy, njx, njy, U, C);

  outputArrayMap(outputFolder, "V", V, 0);
  outputArrayMap(outputFolder, "U", U, 0);
  outputArray(outputFolder, "xc_g", xc_g, 0);
  outputArray(outputFolder, "yc_g", yc_g, 0);
  outputArray(outputFolder, "xc", xc, 0);
  outputArray(outputFolder, "yc", yc, 0);

  outputArray(outputFolder, "Ai", Ai, 0);
  outputArray(outputFolder, "Aj", Aj, 0);
  outputArray(outputFolder, "njx", njx, 0);
  outputArray(outputFolder, "njy", njy, 0);
  outputArray(outputFolder, "nix", nix, 0);
  outputArray(outputFolder, "niy", niy, 0);
  outputArray(outputFolder, "volume", Volume, 0);

  int n = 0;
  cout << "BEFORE" << endl;
  while(time(0,n) < tend && n < 1001)
  {
    for(int k = 0; k < RKORDER; k++)
    {
      setBConU(U_RK, nix, niy, njx, njy, C);
      
      //T = V->Q[pid].cwiseProduct(V->Q[rhoid].cwiseInverse())/R; //computeTemperature
      //setBC(V, nix, niy, njx, njy, T, C);
      //primToCons(U_RK, V);

      computeSourceTerm(S, U_RK, xc, yc, C);

      MUSCL(U_L, U_R, U_B, U_T, U_RK, C);

      compute2dFlux(F, G, U_L, U_R, U_B, U_T, njx, njy, nix, niy, C);

      computeRes(Res, S, F, G, Aj, Ai, Volume, C);

      rungeKutta(U_RK, U, Res, Volume, k, dt, C);

      //consToPrim(V, U_RK);
    }
    U = U_RK;

    dt = computeTimeStepU(Volume, Ai, Aj, nix, niy, njx, njy, U, C);
    //dt = computeTimeStep(Volume, Ai, Aj, nix, niy, njx, njy, V, C);

    time(0,n+1) = time(0,n) + dt;
    n++;

    if(n%C.wint == 0)
    {
      outputArrayMap(outputFolder, "U", U, n);

      outputArrayMap(outputFolder, "F", F, n);
      outputArrayMap(outputFolder, "G", G, n);
      outputArrayMap(outputFolder, "S", S, n);
 
    }

    if(n%C.pint == 0)
      cout << "time  =" << time(0,n) << endl;
  }





 cout << "Made it to End" << endl;

  return 0;
}

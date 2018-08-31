# include "mhdRT.hpp"

int main (int argc, char * argv[])
{
  // Create the structure constants to contain info needed on all processors
  constants C;
  C.f_limiter = 6;
  C.num_ghost = 3;
  C.cfl = 0.5;
  C.nmax = 300;
  C.wint = 50;
  C.pint = 50;
  C.f_mesh = 1;
  C.f_case = 1;
  double dt;
  string output = "output";

  int nx_v, ny_v;
  cout << "Reading In the Mesh Coordinates...." << endl;
  meshSize(&nx_v, &ny_v, C);

  C.nx_c = nx_v-1;
  C.ny_c = ny_v-1;
  double xn[nx_v*ny_v], yn[nx_v*ny_v]; 
  double xc[C.nx_c*C.ny_c], yc[C.nx_c*C.ny_c]; 
  double rc[C.nx_c*C.ny_c], thetc[C.nx_c*C.ny_c]; 
  getCoord(xn, yn, xc, yc, C);
  computeRTheta(rc, thetc, xc, yc, C);

  double xl_g[C.num_ghost*C.ny_c];
  double xr_g[C.num_ghost*C.ny_c];
  double yl_g[C.num_ghost*C.ny_c];
  double yr_g[C.num_ghost*C.ny_c];

  double xb_g[C.num_ghost*C.nx_c];
  double xt_g[C.num_ghost*C.nx_c];
  double yb_g[C.num_ghost*C.nx_c];
  double yt_g[C.num_ghost*C.nx_c];

  cout << "Extrapolating to Ghost Layers...." << endl;
  extrapCopyCoords(xl_g, xr_g, xb_g, xt_g, yl_g, yr_g, yb_g, yt_g, xc, yc, C);

  cout << "Computing Areas & Normal Vectors...." << endl;
  double Aj[nx_v*C.ny_c]; //Number of areas point in the xdir in comp space
  double Ai[ny_v*C.nx_c]; //Number of areas point in the ydir in comp space
  double njx[nx_v*C.ny_c]; double njy[nx_v*C.ny_c]; //Every face has vector
  double nix[ny_v*C.nx_c]; double niy[ny_v*C.nx_c]; //Every face has vector
  computeAreasAndNormalVectors(njx, njy, nix, niy, Aj, Ai, xn, yn, C);

  cout << "Computing Volumes...." << endl;
  double volume[C.nx_c*C.ny_c];
  computeVolume(volume, xn, yn, C);
  outputArray(output, "volume", volume, (sizeof volume/sizeof *volume), 0);
  outputArray(output, "Aj", Aj, (sizeof Aj/sizeof *Aj), 0);
  outputArray(output, "Ai", Ai, (sizeof Ai/sizeof *Ai), 0);
  outputArray(output, "nix", nix, (sizeof nix/sizeof *nix), 0);
  outputArray(output, "niy", niy, (sizeof niy/sizeof *niy), 0);
  outputArray(output, "njx", njx, (sizeof njx/sizeof *njx), 0);
  outputArray(output, "njy", njy, (sizeof njy/sizeof *njy), 0);

  cout << "Initializing...." << endl;
  int nCsize = C.nx_c*C.ny_c*NEQ; // size of NEQ on all cells
  double V[nCsize];  // primitive variables
  double U[nCsize];  // conserved variables
  double Res[nCsize];// Residual 
  double S[nCsize];  // source terms rhou rhov and e wrt g
  initialize(V, rc, thetc, C);
  dt = computeTimeStep(volume, Aj, Ai, njx, njy, nix, niy, V, C);
  //cout << dt << endl;
  prim2Cons(U, V, (sizeof U/sizeof *U));

  double Vl_g[C.num_ghost*C.ny_c*NEQ]; double Vr_g[C.num_ghost*C.ny_c*NEQ];
  double Vt_g[C.num_ghost*C.nx_c*NEQ]; double Vb_g[C.num_ghost*C.nx_c*NEQ];
  double Ul_g[C.num_ghost*C.ny_c*NEQ]; double Ur_g[C.num_ghost*C.ny_c*NEQ];
  double Ut_g[C.num_ghost*C.nx_c*NEQ]; double Ub_g[C.num_ghost*C.nx_c*NEQ];

  cout << "Applying BC...." << endl;
  setBC(Vl_g, Vr_g, Vt_g, Vb_g, njx, njy, nix, niy, V, C);
  prim2Cons(Ul_g, Vl_g, (sizeof Ul_g/ sizeof *Ul_g));
  prim2Cons(Ur_g, Vr_g, (sizeof Ur_g/ sizeof *Ur_g));
  prim2Cons(Ut_g, Vt_g, (sizeof Ut_g/ sizeof *Ut_g));
  prim2Cons(Ub_g, Vb_g, (sizeof Ub_g/ sizeof *Ub_g));

  cout << "MUSCL Extrapolation...." << endl;
  double Ul[nx_v*C.ny_c]; double Ur[nx_v*C.ny_c]; // horiz dir
  double Ub[ny_v*C.nx_c]; double Ut[ny_v*C.nx_c]; // vert dir
  MUSCL(Ul, Ur, Ub, Ut, Ul_g, Ur_g, Ub_g, Ut_g, U, C);

  return 0;

}

/*
//outputArray(output, "V", V, (sizeof V/sizeof *V), 0);
//outputArray(output, "Vl_g", Vl_g, (sizeof Vl_g/sizeof *Vl_g), 0);
//outputArray(output, "U", U, (sizeof U/sizeof *U), 0);
//outputArray(output, "V", V, (sizeof V/sizeof *V), 0);

*/


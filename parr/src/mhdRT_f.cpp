/*
 * Main Function file
 * TS: FEM and GPU
 * Written By: Robert Masti
 * 08/28/2018
 */

# include "mhdRT.hpp"

void computeCoord(int& l, int& r, int& d, int& u, MPI_Comm& com2d)
{
  int rank;
  MPI_Comm_rank(com2d, &rank);
  int coord[2];
  MPI_Cart_coords(com2d, rank, 2, coord);
  MPI_Cart_shift(com2d, 0, 1, &d, &u);
  MPI_Cart_shift(com2d, 1, 1, &l, &r);
}

double computeTimeStepU(
    RowMajorMatrixXd& Volume,      // input - volume of every cell
    RowMajorMatrixXd& Ai,          // input - area in i dir 
    RowMajorMatrixXd& Aj,          // input - area in j dir
    RowMajorMatrixXd& nix,    // input - norm i dir x comp
    RowMajorMatrixXd& niy,    // input - norm i dir y comp
    RowMajorMatrixXd& njx,    // input - norm j dir x comp
    RowMajorMatrixXd& njy,    // input - norm j dir y comp
    Map2Eigen* U,           // input - prim var for speeds
    constants C            // input - constants for the CFL number
    )
{
  int ni = Volume.cols();
  int nj = Volume.rows();

  // need to get the avg normal vectors
  double ni_avg_xhat, ni_avg_yhat;
  double nj_avg_xhat, nj_avg_yhat;
  // get the eigen value which is != MaxSpeed
  double lambda_i, lambda_j;
  double Ai_avg, Aj_avg;

  // add iterables for the cells
  int i_c, j_c;
  int bot = C.num_ghost;// first interior cell in j dir
  int left = C.num_ghost;// first interior cell in i dir
  double dt_min=1e10;
  double maxSpeed;
  double Vmax[NEQ];
  double Umax[NEQ];
  double dt;
  for (int j = 0; j < nj; j++)
  {
    for (int i = 0; i < ni; i++)
    {
      i_c = left + i;
      j_c = bot + j;
      // get the avg norm vecs
      ni_avg_xhat = 0.5*(nix(j,i)+nix(j,i+1));
      ni_avg_yhat = 0.5*(niy(j,i)+niy(j,i+1));
      nj_avg_xhat = 0.5*(njx(j,i)+njx(j+1,i));
      nj_avg_yhat = 0.5*(njy(j,i)+njy(j+1,i));

      // get avg area
      Ai_avg = 0.5*(Ai(j,i) + Ai(j,i+1));
      Aj_avg = 0.5*(Aj(j,i) + Aj(j+1,i));

      // compute max eigen values

      for(int eq = 0; eq < NEQ; eq++)
        Umax[eq] = U->Q[eq](j_c,i_c);

      cons2Prim(Vmax, Umax, NEQ);
      maxSpeed = computeMaxSpeed(Vmax,0);
      lambda_i = abs(Vmax[uid]*ni_avg_xhat 
          + Vmax[vid]*ni_avg_yhat) + maxSpeed;

      maxSpeed = computeMaxSpeed(Vmax,1);
      lambda_j = abs(Vmax[uid]*nj_avg_xhat 
          + Vmax[vid]*nj_avg_yhat) + maxSpeed;

      // compute time step non uniform grid
      dt = C.cfl*Volume(j,i) / (lambda_i*Ai_avg + lambda_j*Aj_avg);
      if (dt < dt_min)
        dt_min = dt;
    }
  }
  return dt_min;
}



void setBConU(Map2Eigen* U, const RowMajorMatrixXd& nix, const RowMajorMatrixXd& niy, const RowMajorMatrixXd& njx, const RowMajorMatrixXd njy, constants C)
{
  int ic, jc;
  int nj = U->Q[rhoid].rows();
  int ni = U->Q[rhoid].cols();
  int igr;
  int igl;
  int icr;
  int icl;

  // left and right use periodic BC
  icl = C.num_ghost;
  icr = ni-C.num_ghost-1;
  for(int i = 0; i < C.num_ghost; i++)
  {

    igr = (ni-C.num_ghost)+i;
    igl = (C.num_ghost-1)-i;

    // Left g copy right interior
    //indepndent
    U->Q[rhoid].col(igl) = U->Q[rhoid].col(icr-i);
    //cout << U->Q[rhoid].col(icr-i) << endl;
    U->Q[wid].col(igl) = U->Q[wid].col(icr-i);
    U->Q[pid].col(igl) = U->Q[pid].col(icr-i);
    U->Q[bzid].col(igl) = U->Q[bzid].col(icr-i);

    //direction flip
    U->Q[uid].col(igl) = -1*U->Q[uid].col(icr-i);
    U->Q[vid].col(igl) = -1*U->Q[vid].col(icr-i);
    U->Q[bxid].col(igl) = -1*U->Q[bxid].col(icr-i);
    U->Q[byid].col(igl) = -1*U->Q[byid].col(icr-i);

    // Right g copy left interior
    //indepndent
    U->Q[rhoid].col(igr) = U->Q[rhoid].col(icl+i);
    U->Q[wid].col(igr) = U->Q[wid].col(icl+i);
    U->Q[pid].col(igr) = U->Q[pid].col(icl+i);
    U->Q[bzid].col(igr) = U->Q[bzid].col(icl+i);

    //direction flip
    U->Q[uid].col(igr) = -1*U->Q[uid].col(icl+i);
    U->Q[vid].col(igr) = -1*U->Q[vid].col(icl+i);
    U->Q[bxid].col(igr) = -1*U->Q[bxid].col(icl+i);
    U->Q[byid].col(igr) = -1*U->Q[byid].col(icl+i);
  }

  // top bottom slipwall
  int jfb, jft;
  jfb = 0; jft = njx.rows()-1; //top face
  int jgb, jgt;
  int iff;

  int jcb, jct; 
  double uvel, vvel, bx, by, nx,ny;
  for(int j = 0; j < C.num_ghost; j++)
  {
    for(int i = C.num_ghost; i < ni-C.num_ghost; i++)
    {

      iff = i-C.num_ghost;
      //for(int i = C.num_ghost; i < ni-C.num_ghost; i++)
      jcb = C.num_ghost + j;
      jct = nj-C.num_ghost-1 - j;

      jgt = (nj-C.num_ghost)+j;
      jgb = (C.num_ghost-1)-j;  
      //bottom 

      // Get vectors

      nx = njx(jfb,iff); 
      ny = njy(jfb,iff);
      uvel = U->Q[uid](jcb+j,i);
      vvel = U->Q[vid](jcb+j,i);
      bx = U->Q[bxid](jcb+j,i);
      by = U->Q[byid](jcb+j,i);

      //horizontal boundary lower
      //no penetrate
      U->Q[uid](jgb,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx); 
      U->Q[vid](jgb,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
      U->Q[bxid](jgb,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx); 
      U->Q[byid](jgb,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
      // extrapolate
      U->Q[rhoid](jgb,i) = 2.0*U->Q[rhoid](jgb+1,i) - U->Q[rhoid](jgb+2,i);
      U->Q[wid](jgb,i) = 2.0*U->Q[wid](jgb+1,i) - U->Q[wid](jgb+2,i);
      U->Q[pid](jgb,i) = 2.0*U->Q[pid](jgb+1,i) - U->Q[pid](jgb+2,i);
      U->Q[bzid](jgb,i) = 2.0*U->Q[bzid](jgb+1,i) - U->Q[bzid](jgb+2,i);

      /*
         U->Q[rhoid](jgb,i) = U->Q[rhoid](jgb+1,i);// - U->Q[rhoid](jgb+2,i);
         U->Q[wid](jgb,i) = U->Q[wid](jgb+1,i);// - U->Q[wid](jgb+2,i);
         U->Q[pid](jgb,i) = U->Q[pid](jgb+1,i);// - U->Q[pid](jgb+2,i);
         U->Q[bzid](jgb,i) = U->Q[bzid](jgb+1,i);// - U->Q[bzid](jgb+2,i);

*/

      // horizontal boundary upper
      nx = njx(jft,iff); 
      ny = njy(jft,iff);
      uvel = U->Q[uid](jct-j,i);
      vvel = U->Q[vid](jct-j,i);
      bx = U->Q[bxid](jct-j,i);
      by = U->Q[byid](jct-j,i);

      //no penetrate
      U->Q[uid](jgt,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx); 
      U->Q[vid](jgt,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
      U->Q[bxid](jgt,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx); 
      U->Q[byid](jgt,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);

      // extrapolate
      U->Q[rhoid](jgt,i) = 2.0*U->Q[rhoid](jgt-1,i) - U->Q[rhoid](jgt-2,i);
      U->Q[wid](jgt,i) = 2.0*U->Q[wid](jgt-1,i) - U->Q[wid](jgt-2,i);
      U->Q[pid](jgt,i) = 2.0*U->Q[pid](jgt-1,i) - U->Q[pid](jgt-2,i);
      U->Q[bzid](jgt,i) = 2.0*U->Q[bzid](jgt-1,i) - U->Q[bzid](jgt-2,i);

      /*
         U->Q[rhoid](jgt,i) = U->Q[rhoid](jgt-1,i);// - U->Q[rhoid](jgt-2,i);
         U->Q[wid](jgt,i) = U->Q[wid](jgt-1,i);// - U->Q[wid](jgt-2,i);
         U->Q[pid](jgt,i) = U->Q[pid](jgt-1,i);// - U->Q[pid](jgt-2,i);
         U->Q[bzid](jgt,i) =U->Q[bzid](jgt-1,i);// - U->Q[bzid](jgt-2,i);

*/

    }
  }
}

void rungeKutta(
    // This function performs a runge kutta updating scheme based on the k
    Map2Eigen* U_RK,        // output - Consvar with RK dt stepping  
    Map2Eigen* U,           // input - Consvar before
    Map2Eigen* Res,         // input - Residual 
    RowMajorMatrixXd& Volume,      // input - Volume matrix
    int k,                 // input - iteration
    double dt,         // input - dt minimum
    constants C            // input - constants 
    ) 
{
  int ni = Res->Q[rhoid].cols(); // number int cells i dir
  int nj = Res->Q[rhoid].rows(); // number in cells j dir
  int bot = C.num_ghost;
  int left = C.num_ghost;
  if (Volume.rows() != nj || Volume.cols() != ni)
  {
    cerr << "ERROR: Dimension Mismatch in RK?!?!" << endl;
    exit(1);
  } 
  double a[4];
  if (RKORDER == 1) // 1 stage RK 
    a[0] = 1.0; a[1]=a[2]=a[3]=0;
  if (RKORDER == 2) // 2 stage RK
  {
    a[0] = 0.5; a[1] = 1.0; a[2]=a[3]=0;
  }
  if (RKORDER == 4) // 4 stage RK
  {
    a[0] = 0.25; a[1] = 1.0/3.0; a[2] = 0.5; a[3] = 1.0;
  }
  if (RKORDER == 3 || RKORDER > 4) // no sense in doing 3
  {
    cerr << "ERROR: RK Order Not Available?!?!" << endl;
    exit(1);
  } 
  // set iteration for big matrix
  int i_c, j_c;
  for (int j = 0; j < nj; j++)
    for (int i = 0; i < ni; i++)
    {
      i_c = left + i;
      j_c = bot + j;
      for (int eq = 0; eq < NEQ; eq++) 
        U_RK->Q[eq](j_c,i_c) = U->Q[eq](j_c,i_c) - a[k]*dt*Res->Q[eq](j,i)/Volume(j,i);
    }
}

void computeRes(Map2Eigen* Res, const Map2Eigen* S, const Map2Eigen* F, const Map2Eigen* G, const RowMajorMatrixXd& Aj, const RowMajorMatrixXd& Ai, const RowMajorMatrixXd& Volume, constants C)
{
  for (int j = 0; j < Res->Q[rhoid].rows(); j++)
    for (int i = 0; i < Res->Q[rhoid].cols(); i++)
      // res = right + top - left - bottom - source
      for (int eq = 0; eq < NEQ; eq++)
        Res->Q[eq](j,i) = F->Q[eq](j,i+1)*Ai(j,i+1) 
          + G->Q[eq](j+1,i)*Aj(j+1,i)
          - F->Q[eq](j,i)*Ai(j,i) 
          - G->Q[eq](j,i)*Aj(j,i) 
          - S->Q[eq](j,i)*Volume(j  ,i);
}

void computeSourceTerm(Map2Eigen* S, Map2Eigen* U, const RowMajorMatrixXd& xc, const RowMajorMatrixXd& yc, constants C)
{
  int nj = xc.rows();
  int ni = xc.cols();
  int ig, jg;

  /*
     S->Q[rhoid] = RowMajorMatrixXd::Constant(nj,ni,0.0); // sets all values to rho 
     S->Q[wid] = RowMajorMatrixXd::Constant(nj,ni,0.0); // sets all values to rho 
     S->Q[bxid] = RowMajorMatrixXd::Constant(nj,ni,0.0); // sets all values to rho 
     S->Q[byid] = RowMajorMatrixXd::Constant(nj,ni,0.0); // sets all values to rho 
     S->Q[bzid] = RowMajorMatrixXd::Constant(nj,ni,0.0); // sets all values to rho 
     */
  double thet; 
  double g =ACCEL;// 0.1;
  for(int j = 0; j < nj; j++)
  {
    jg = j + C.num_ghost;
    for(int i = 0; i < ni; i++)
    {
      ig = i + C.num_ghost;
      thet = atan2(yc(j,i), xc(j,i));
      S->Q[uid](j,i) = -1*(U->Q[rhoid](jg,ig))*(g)*cos(thet);
      S->Q[vid](j,i) = -1*(U->Q[rhoid](jg,ig))*(g)*sin(thet);
      S->Q[pid](j,i) = ((S->Q[uid](j,i)*-1*g*cos(thet)) + (S->Q[vid](j,i)*-1*g*sin(thet)));//U->Q[rhoid](jg,ig) ;
    }
  }
}

void compute2dFlux(Map2Eigen* F, Map2Eigen* G,  Map2Eigen* U_L ,  Map2Eigen* U_R,  Map2Eigen* U_B,  Map2Eigen* U_T, RowMajorMatrixXd& njx, RowMajorMatrixXd& njy, RowMajorMatrixXd& nix, RowMajorMatrixXd& niy, constants C)
{
  int nj = F->Q[rhoid].rows();
  int ni = G->Q[rhoid].cols();
  double ul[NEQ], ur[NEQ], ub[NEQ], ut[NEQ];
  double FFLUX[NEQ];
  double GFLUX[NEQ];

  for(int j = 0; j < nj; j++)
  {
    for(int i = 0; i < ni; i++)
    {
      for(int eq = 0; eq < NEQ; eq++)
      {
        ul[eq] = U_L->Q[eq](j,i);
        ur[eq] = U_R->Q[eq](j,i);
        ub[eq] = U_B->Q[eq](j,i);
        ut[eq] = U_T->Q[eq](j,i);
      }

      computeFluxHLL(FFLUX, ul, ur, nix(j,i), niy(j,i), 0);
      computeFluxHLL(GFLUX, ub, ut, njx(j,i), njy(j,i), 1);

      //computeFluxVL(FFLUX, ul, ur, nix(j,i), niy(j,i));
      //computeFluxVL(GFLUX, ub, ut, njx(j,i), njy(j,i));

      //computeFluxRoe(FFLUX, ul, ur, nix(j,i), niy(j,i));
      //computeFluxRoe(GFLUX, ub, ut, njx(j,i), njy(j,i));

      //computeFluxHLLD(FFLUX, ul, ur, nix(j,i), niy(j,i), 0);
      //computeFluxHLLD(GFLUX, ub, ut, njx(j,i), njy(j,i), 1);



      for(int eq = 0; eq < NEQ; eq++)
      {
        F->Q[eq](j,i) = FFLUX[eq];
        G->Q[eq](j,i) = GFLUX[eq];
      }
    }
  }

  // fix right wall
  for(int j = 0; j < nj; j++)
  {
    for(int eq = 0; eq < NEQ; eq++)
    {
      ul[eq] = U_L->Q[eq](j,ni);
      ur[eq] = U_R->Q[eq](j,ni);

    }

    computeFluxHLL(FFLUX, ul, ur, nix(j,ni), niy(j,ni), 0);
    //computeFluxVL(FFLUX, ul, ur, nix(j,ni), niy(j,ni));
    //computeFluxRoe(FFLUX, ul, ur, nix(j,ni), niy(j,ni));
    //computeFluxHLLD(FFLUX, ul, ur, nix(j,ni), niy(j,ni), 0);
    for(int eq = 0; eq < NEQ; eq++)
      F->Q[eq](j,ni) = FFLUX[eq];
  }
  // fix left wall

  for(int i = 0; i < ni; i++)
  {
    for(int eq = 0; eq < NEQ; eq++)
    {
      ub[eq] = U_B->Q[eq](nj,i);
      ut[eq] = U_T->Q[eq](nj,i);

    }
    computeFluxHLL(GFLUX, ub, ut, njx(nj,i), njy(nj,i), 1);
    //computeFluxVL(GFLUX, ub, ut, njx(nj,i), njy(nj,i));
    //computeFluxRoe(GFLUX, ub, ut, njx(nj,i), njy(nj,i));
    //computeFluxHLLD(GFLUX, ub, ut, njx(nj,i), njy(nj,i), 1);
    for(int eq = 0; eq < NEQ; eq++)
      G->Q[eq](nj,i) = GFLUX[eq];
  }
}

void computeFluxHLL(double F[], double UA[], double UB[], double& nxhat, double& nyhat, int ForG)
{



  double CA = computeMaxSpeed(UA, ForG);
  double CB = computeMaxSpeed(UB, ForG);

  double u_A, u_B;

  u_A = (UA[uid]*nxhat + UA[vid]*nyhat)/UA[rhoid]; 
  u_B = (UB[uid]*nxhat + UB[vid]*nyhat)/UB[rhoid]; 

  double lambdaA, lambdaB;
  lambdaA = mymin(u_A, u_B) - mymax(CA, CB);

  lambdaB = mymax(u_A, u_B) + mymax(CA, CB);

  double FA[NEQ], FB[NEQ];

  if (ForG == 0)
  {
    fFlux(FA, UA, nxhat, nyhat);
    fFlux(FB, UB, nxhat, nyhat);
  }
  else if (ForG == 1)
  {
    gFlux(FA, UA, nxhat, nyhat);
    gFlux(FB, UB, nxhat, nyhat);
  }
  else
  {
    cerr << "ERROR: Wrong index for direction" << endl;
    exit(-1);
  }

  //for(int eq=0; eq < NEQ; eq++)
  if (lambdaA > 0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = FA[eq];
  else if (lambdaB < 0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = FB[eq];
  else if (lambdaB >=0 && lambdaA <=0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = (lambdaB*FA[eq] - lambdaA*FB[eq]+lambdaA*lambdaB*(UB[eq]-UA[eq]))/(lambdaB-lambdaA);
  else
  {
    cerr << "ERROR: HLL Flux Calculation Lambda Unphysical" << endl;
    exit(-1);
  }
}
void fFlux(double F[], double U[], double nxhat, double nyhat)
{
  double e;
  double V[NEQ];
  cons2Prim(V, U, NEQ);


  double uHat = V[uid]*nxhat+V[vid]*nyhat;
  double bxHat = V[bxid]*nxhat+V[byid]*nyhat;

  F[rhoid] = V[rhoid]*uHat; 
  F[uid] = V[rhoid]*V[uid]*uHat - V[bxid]*bxHat/MU0 + V[pid]*nxhat + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])*nxhat/(2*MU0);

  F[vid] = V[rhoid]*uHat*V[vid]- bxHat*V[byid]/MU0+ V[pid]*nyhat + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])*nyhat/(2*MU0);

  F[wid] = V[rhoid]*uHat*V[wid] - bxHat*V[bzid]/MU0;

  e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

  F[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*uHat - (1/MU0)*bxHat*(V[uid]*V[bxid]+V[vid]*V[byid]+V[wid]*V[bzid]);

  F[bxid] = uHat*V[bxid]-V[uid]*bxHat;

  F[byid] = uHat*V[byid]-V[vid]*bxHat;

  F[bzid] = uHat*V[bzid]-V[wid]*bxHat;

  //for(int eq=0; eq < NEQ; eq++)
  //F[eq] = 0.0;
  /*
     F[rhoid] = V[rhoid]*V[uid]; 
     F[uid] = V[rhoid]*V[uid]*V[uid] - (V[bxid]*V[bxid])/MU0 + V[pid] + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])/(2*MU0);

     F[vid] = V[rhoid]*V[uid]*V[vid]- V[bxid]*V[byid]/MU0;

     F[wid] = V[rhoid]*V[uid]*V[wid] - V[bxid]*V[bzid]/MU0;

     e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

     F[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*V[uid] - (1/MU0)*V[bxid]*(V[uid]*V[bxid]+V[vid]*V[byid]+V[wid]*V[bzid]);

     F[bxid] = 0;

     F[byid] = V[uid]*V[byid]-V[vid]*V[bxid];

     F[bzid] = V[uid]*V[bzid]-V[wid]*V[bxid];
     */
}
void gFlux(double G[], double U[], double nxhat, double nyhat)
{
  double e;
  double V[NEQ];
  cons2Prim(V, U, NEQ);

  double vHat = V[uid]*nxhat+V[vid]*nyhat;
  double byHat = V[bxid]*nxhat+V[byid]*nyhat;


  G[rhoid] = V[rhoid]*vHat;

  G[uid] = V[rhoid]*V[uid]*vHat- V[bxid]*byHat/MU0 + V[pid]*nxhat + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])*nxhat/(2*MU0);

  G[vid] = V[rhoid]*V[vid]*vHat - V[byid]*byHat/MU0 + V[pid]*nyhat + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])*nyhat/(2*MU0);

  G[wid] = V[rhoid]*vHat*V[wid]-byHat*V[bzid]/MU0;

  e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

  G[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*vHat - (1/MU0)*byHat*(V[uid]*V[bxid]+V[vid]*V[byid]+ V[wid]*V[bzid]);

  G[bxid] = vHat*V[bxid]-V[uid]*byHat;

  G[byid] = vHat*V[byid]-V[vid]*byHat;

  G[bzid] = vHat*V[bzid]-V[wid]*byHat;

  //for(int eq=0; eq < NEQ; eq++)
  //G[eq] = 0.0;
  /*
     G[rhoid] = V[rhoid]*V[vid];

     G[uid] = V[rhoid]*V[uid]*V[vid]- V[bxid]*V[byid]/MU0;

     G[vid] = V[rhoid]*V[vid]*V[vid] - (V[byid]*V[byid])/MU0 + V[pid] + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])/(2*MU0);

     G[wid] = V[rhoid]*V[vid]*V[wid]-V[byid]*V[bzid]/MU0;

     e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

     G[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*V[vid] - (1/MU0)*V[byid]*(V[uid]*V[bxid]+V[vid]*V[byid]+ V[wid]*V[bzid]);

     G[bxid] = V[vid]*V[bxid]-V[uid]*V[byid];

     G[byid] = 0;

     G[bzid] = V[vid]*V[bzid]-V[wid]*V[byid];
     */
}




void MUSCL(Map2Eigen* U_L, Map2Eigen* U_R, Map2Eigen* U_B, Map2Eigen* U_T,const Map2Eigen* U, constants C)
{


  int ni = U_B->Q[rhoid].cols(); //ncx
  int nj = U_L->Q[rhoid].rows();//ncy

  int ig, jg;
  double theta_L, theta_B, theta_R, theta_T;
  double denom_L, denom_B, denom_R, denom_T;
  double r_L, r_B, r_R, r_T;
  double delta = 1.0e-6;

  /*
     printf("\n nj = %d, ni = %d \n", nj, ni);
     printf("\n nj = %ld, ni = %ld \n", U_B->Q[rhoid].rows() , U_B->Q[rhoid].cols());
     printf("\n nj = %ld, ni = %ld \n", U_R->Q[rhoid].rows(), U_R->Q[rhoid].cols());
     printf("\n nj = %ld, ni = %ld \n", U_L->Q[rhoid].rows(), U_L->Q[rhoid].cols());
     printf("\n nj = %ld, ni = %ld \n", U_T->Q[rhoid].rows(), U_T->Q[rhoid].cols());
     printf("\n nj = %ld, ni = %ld \n", U->Q[rhoid].rows(), U->Q[rhoid].cols());
     */
  for(int j = 0; j < nj; j++)
  {
    for(int i = 0; i < ni; i++)
    {
      ig = i+C.num_ghost;
      jg = j+C.num_ghost;

      for(int eq=0; eq < NEQ; eq++)
      {
        // left right

        denom_L = U->Q[eq](jg, ig) - U->Q[eq](jg, ig-1);
        denom_R = U->Q[eq](jg, ig+1) - U->Q[eq](jg, ig-1);
        denom_L = SIGN(1,denom_L)*mymax(delta, abs(denom_L));    
        denom_R = SIGN(1,denom_R)*mymax(delta, abs(denom_R));    
        r_L = (U->Q[eq](jg,ig-1) - U->Q[eq](jg,ig-2))/denom_L;
        r_R = (U->Q[eq](jg,ig) - U->Q[eq](jg,ig-1))/denom_R;
        theta_L = limiter(r_L, C);
        theta_R = limiter(r_R, C);

        denom_B = U->Q[eq](jg, ig) - U->Q[eq](jg-1, ig);
        denom_T = U->Q[eq](jg+1, ig) - U->Q[eq](jg-1, ig);
        denom_B = SIGN(1,denom_B)*mymax(delta, abs(denom_B));    
        denom_T = SIGN(1,denom_T)*mymax(delta, abs(denom_T));    
        r_B = (U->Q[eq](jg-1,ig) - U->Q[eq](jg-2,ig))/denom_B;
        r_T = (U->Q[eq](jg,ig) - U->Q[eq](jg-1,ig))/denom_T;
        theta_B = limiter(r_B, C);
        theta_T = limiter(r_T, C);

        U_R->Q[eq](j,i) = (U->Q[eq](jg,ig))-0.5*theta_R*((U->Q[eq](jg,ig+1))-(U->Q[eq](jg,ig)));             
        //(U_R->Q_raw[(j+i)*NEQ+eq]) = 100.0;

        U_L->Q[eq](j,i) =  (U->Q[eq](jg,ig-1))+0.5*theta_L*((U->Q[eq](jg,ig))-(U->Q[eq](jg,ig-1)));

        U_T->Q[eq](j,i) = (U->Q[eq](jg,ig))+0.5*theta_T*((U->Q[eq](jg+1,ig))-(U->Q[eq](jg,ig)));

        U_B->Q[eq](j,i) =  (U->Q[eq](jg-1,ig))+0.5*theta_B*((U->Q[eq](jg,ig))-(U->Q[eq](jg-1,ig)));

      }
    }
  }

  // add the last layer left right
  //
  for(int j = 0; j < nj; j++)
  {
    ig = ni+C.num_ghost;
    jg = j+C.num_ghost;
    for(int eq = 0; eq < NEQ; eq++)
    {
      denom_L = U->Q[eq](jg, ig) - U->Q[eq](jg, ig-1);
      denom_R = U->Q[eq](jg, ig+1) - U->Q[eq](jg, ig-1);
      denom_L = SIGN(1,denom_L)*mymax(delta, abs(denom_L));    
      denom_R = SIGN(1,denom_R)*mymax(delta, abs(denom_R));    
      r_L = (U->Q[eq](jg,ig-1) - U->Q[eq](jg,ig-2))/denom_L;
      r_R = (U->Q[eq](jg,ig) - U->Q[eq](jg,ig-1))/denom_R;
      theta_L = limiter(r_L, C);
      theta_R = limiter(r_R, C);

      U_L->Q[eq](j, ni) =  U->Q[eq](jg,ig-1)+0.5*theta_L*(U->Q[eq](jg,ig)-U->Q[eq](jg,ig-1));
      U_R->Q[eq](j, ni) = U->Q[eq](jg,ig)-0.5*theta_R*(U->Q[eq](jg,ig+1)-U->Q[eq](jg,ig));             
    }
  }
  // add the last layer top bottom
  for(int i = 0; i < ni; i++)
  {
    ig = i+C.num_ghost;
    jg = nj+C.num_ghost;
    for(int eq = 0; eq < NEQ; eq++)
    {
      denom_B = U->Q[eq](jg, ig) - U->Q[eq](jg-1, ig);
      denom_T = U->Q[eq](jg+1, ig) - U->Q[eq](jg-1, ig);
      denom_B = SIGN(1,denom_B)*mymax(delta, abs(denom_B));    
      denom_T = SIGN(1,denom_T)*mymax(delta, abs(denom_T));    
      r_B = (U->Q[eq](jg-1,ig) - U->Q[eq](jg-2,ig))/denom_B;
      r_T = (U->Q[eq](jg,ig) - U->Q[eq](jg-1,ig))/denom_T;
      theta_B = limiter(r_B, C);
      theta_T = limiter(r_T, C);

      U_B->Q[eq](nj, i) =  U->Q[eq](jg-1,ig)+0.5*theta_B*(U->Q[eq](jg,ig)-U->Q[eq](jg-1,ig));
      U_T->Q[eq](nj,i) = U->Q[eq](jg,ig)-0.5*theta_T*(U->Q[eq](jg+1,ig)-U->Q[eq](jg,ig));             
    }
  }
}


double limiter(double& r, constants C)
{
  double theta;
  switch (C.f_limiter)
  {
    case 1:
      theta = 1.0;
    case 2:
      theta = (r+abs(r))/(1+abs(r));
    case 3:
      theta = (r+r*r)/(1+r*r);
    case 4:
      theta = 1.5*(r*r + r) / (1 + r + r*r);
    case 5:
      theta = mymax(0,mymin(2.0*r,mymin(0.5*(1+r),2)));
    case 6:
      theta = mymax(0,mymin(1.0,r));
    case 7:
      theta = mymax(0,mymax(mymin(2.0*r,1),mymin(r,2)));
  }
  return theta;
}            


double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}


double computeTimeStep(
    RowMajorMatrixXd& Volume,      // input - volume of every cell
    RowMajorMatrixXd& Ai,          // input - area in i dir 
    RowMajorMatrixXd& Aj,          // input - area in j dir
    RowMajorMatrixXd& nix,    // input - norm i dir x comp
    RowMajorMatrixXd& niy,    // input - norm i dir y comp
    RowMajorMatrixXd& njx,    // input - norm j dir x comp
    RowMajorMatrixXd& njy,    // input - norm j dir y comp
    Map2Eigen* V,           // input - prim var for speeds
    constants C            // input - constants for the CFL number
    )
{

  int ni = Volume.cols();
  int nj = Volume.rows();

  // need to get the avg normal vectors
  double ni_avg_xhat, ni_avg_yhat;
  double nj_avg_xhat, nj_avg_yhat;
  // get the eigen value which is != MaxSpeed
  double lambda_i, lambda_j;
  double Ai_avg, Aj_avg;

  // add iterables for the cells
  int i_c, j_c;
  int bot = C.num_ghost;// first interior cell in j dir
  int left = C.num_ghost;// first interior cell in i dir
  double dt_min=1e10;
  double maxSpeed;
  double Vmax[NEQ];
  double dt;
  for (int j = 0; j < nj; j++)
  {
    for (int i = 0; i < ni; i++)
    {
      i_c = left + i;
      j_c = bot + j;
      // get the avg norm vecs
      ni_avg_xhat = 0.5*(nix(j,i)+nix(j,i+1));
      ni_avg_yhat = 0.5*(niy(j,i)+niy(j,i+1));
      nj_avg_xhat = 0.5*(njx(j,i)+njx(j+1,i));
      nj_avg_yhat = 0.5*(njy(j,i)+njy(j+1,i));

      // get avg area
      Ai_avg = 0.5*(Ai(j,i) + Ai(j,i+1));
      Aj_avg = 0.5*(Aj(j,i) + Aj(j+1,i));

      // compute max eigen values

      for(int eq = 0; eq < NEQ; eq++)
        Vmax[eq] = V->Q[eq](j_c,i_c);

      maxSpeed = computeMaxSpeed(Vmax, 0);
      lambda_i = abs(V->Q[uid](j_c,i_c)*ni_avg_xhat 
          + V->Q[vid](j_c,i_c)*ni_avg_yhat) + maxSpeed;
      maxSpeed = computeMaxSpeed(Vmax, 1);
      lambda_j = abs(V->Q[uid](j_c,i_c)*nj_avg_xhat 
          + V->Q[vid](j_c,i_c)*nj_avg_yhat) + maxSpeed;

      // compute time step non uniform grid
      dt = C.cfl*Volume(j,i) / (lambda_i*Ai_avg + lambda_j*Aj_avg);
      if (dt < dt_min)
        dt_min = dt;
    }
  }
  return dt_min;
}

double computeMaxSpeed(double V[], int ForG)                                         
{
  double rho, u, v, w, p, bx, by, bz;
  rho = V[rhoid];
  u = V[uid];
  v = V[vid];
  w = V[wid];
  p = V[pid];
  bx = V[bxid];
  by = V[byid];
  bz = V[bzid];
  double cs, ca, cadir;
  double bmag, vmag;
  double maxSpeed;
  bmag = sqrt(bx*bx+by*by+bz*bz);
  vmag = sqrt(u*u+v*v+w*w);
  cs = sqrt(p*GAMMA/rho);
  ca = bmag/sqrt(MU0*rho);


  if (ForG == 0)
    cadir = bx*bx/(MU0*rho);
  else
    cadir = by*by/(MU0*rho);


  maxSpeed = sqrt(0.5*(ca*ca+cs*cs+sqrt((ca*ca+cs*cs)*(ca*ca+cs*cs)-4*cs*cs*cadir*cadir)));// fast magnetosonic speed

  return maxSpeed;
}                                                                                       



void slipwallBC(
    // This function will allow the normal velocity from the boundary to reflect but simultan
    // eously keep the parrallel to the boundary vel unchanged
    Map2Eigen* V,              // output - Prim var  
    const int Begin[],         // input - Beginning index coord
    const int End[],           // input - Ending index coord 
    const RowMajorMatrixXd& nix,  // input - norm vec i dir x comp
    const RowMajorMatrixXd& niy,  // input - norm vec i dir y comp
    const RowMajorMatrixXd& njx,  // input - norm vec j dir x comp
    const RowMajorMatrixXd& njy,  // input - norm vec j dir y comp
    RowMajorMatrixXd& T,         // input - temperature at cells
    constants C                // input - constants for num ghosts
    )
{
  int ni_g = V->Q[rhoid].cols();
  int nj_g = V->Q[rhoid].rows();

  int ni = ni_g - 2*C.num_ghost;
  int nj = nj_g - 2*C.num_ghost;

  int i_in_Begin = Begin[1] + C.num_ghost; // 1st interior i dir
  int j_in_Begin = Begin[0] + C.num_ghost; // 1st interior j dir
  int i_in_End = End[1] + C.num_ghost; // shifted by numghost
  int j_in_End = End[0] + C.num_ghost;
  int I, J;
  int sign;
  int BC_index;
  double nx, ny, uvel, vvel;
  double Temp;

  if (i_in_Begin == i_in_End) // Then we are along a vertical boundary loop over j
  {
    if (i_in_Begin == C.num_ghost) // detect whether the ghost cells are on the left (i=3)
    {
      sign = -1; // the face points in the - direction if on the left
      BC_index = Begin[1];
    }
    else if (i_in_Begin == ni + C.num_ghost - 1) // now its on the right wall so sign + and 
    {
      sign = 1;
      BC_index = Begin[1] + 1;
    }
    else
    {
      cerr << "ERROR: Not Applying the Correct SlipWall BC's col!!" << endl;
      exit(1);
    }

    // LOOP OVER THE J-Dir for a vertical boundary!!!!!!!!!!!!!
    for (int j = 0; j <= End[0] - Begin[0]; j++)
    {
      for (int i = 0; i < C.num_ghost; i++)
      {
        I = i_in_Begin + sign*1; // gives the first coln inside the ghost cell
        J = j_in_Begin; // this give first j index
        nx = nix(Begin[0]+j, BC_index); // get the normal vec comp in the xdir
        ny = niy(Begin[0]+j, BC_index); // get normal vec comp in the y dir
        uvel = V->Q[uid](J+j, i_in_Begin - sign*i); // sign says okay push left for g
        vvel = V->Q[vid](J+j, i_in_Begin - sign*i); // or push right for g
        Temp = T(J+j, i_in_Begin - sign*i); // temperature 

        // get the velocity into comp coords
        V->Q[uid](J+j, I+sign*i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        V->Q[vid](J+j, I+sign*i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);

        // extrapolate P
        V->Q[pid](J+j, I+sign*i) =  2.0*V->Q[pid](J+j, I+sign*(i-1))-V->Q[pid](J+j, I+sign*(i-2));
        V->Q[bxid](J+j, I+sign*i) =  2.0*V->Q[bxid](J+j, I+sign*(i-1))-V->Q[bxid](J+j, I+sign*(i-2));
        V->Q[byid](J+j, I+sign*i) =  2.0*V->Q[byid](J+j, I+sign*(i-1))-V->Q[byid](J+j, I+sign*(i-2));
        V->Q[bzid](J+j, I+sign*i) =  2.0*V->Q[bzid](J+j, I+sign*(i-1))-V->Q[bzid](J+j, I+sign*(i-2));
        V->Q[wid](J+j, I+sign*i) =  2.0*V->Q[wid](J+j, I+sign*(i-1))-V->Q[wid](J+j, I+sign*(i-2));
        T(J+j, I+sign*i) = Temp; // copy temperature extrapolate pressure then compute rho
        V->Q[rhoid](J+j, I+sign*i) = V->Q[pid](J+j, I+sign*i) / (R*T(J+j, I+sign*i));
      }
    }
  }
  else if (j_in_Begin == j_in_End) // now have horizontal boundary
  {
    if (j_in_Begin == C.num_ghost) // at lowere
    {
      sign = -1;
      BC_index = Begin[0];
    }
    else if ( j_in_Begin == nj + C.num_ghost-1 ) // top boundary
    {
      sign = 1;
      BC_index = Begin[0] + 1;
    }
    else
    {
      cerr << "ERROR: Not Applying the Correct SlipWall BC's row!!" << endl;
      exit(1);
    }
    for (int j = 0 ; j < C.num_ghost; j++)
    {
      for (int i = 0 ; i <= End[1]-Begin[1]; i++)
      {
        I = i_in_Begin;          // I is the looping dir
        J = j_in_Begin + sign*1; // j is the num ghost dir
        nx = njx(BC_index, Begin[1]+i);
        ny = njy(BC_index, Begin[1]+i);
        uvel = V->Q[uid](j_in_Begin - sign*j, I+i); // copy val
        vvel = V->Q[vid](j_in_Begin - sign*j, I+i); // copy val
        Temp = T(j_in_Begin - sign*j, I+i); // copy val
        V->Q[uid](J+sign*j, I+i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        V->Q[vid](J+sign*j, I+i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        V->Q[pid](J+sign*j, I+i) = 2.0*V->Q[pid](J+sign*(j-1), I+i) - V->Q[pid](J+sign*(j-2), I+i);
        V->Q[wid](J+sign*j, I+i) = 2.0*V->Q[wid](J+sign*(j-1), I+i) - V->Q[wid](J+sign*(j-2), I+i);
        V->Q[bxid](J+sign*j, I+i) = 2.0*V->Q[bxid](J+sign*(j-1), I+i) - V->Q[bxid](J+sign*(j-2), I+i);
        V->Q[byid](J+sign*j, I+i) = 2.0*V->Q[byid](J+sign*(j-1), I+i) - V->Q[byid](J+sign*(j-2), I+i);
        V->Q[bzid](J+sign*j, I+i) = 2.0*V->Q[bzid](J+sign*(j-1), I+i) - V->Q[bzid](J+sign*(j-2), I+i);
        T(J+sign*j, I+i) = Temp;
        V->Q[rhoid](J+sign*j, I+i) = V->Q[pid](J+sign*j, I+i) / (R*T(J+sign*j, I+i));
      }
    }
  }
  else
  {
    cerr << "ERROR: Index Incorrect for SlipWall !!" << endl;
    exit(1);
  }
}

void setBC(
    // Apply BC for cases 2 and 3, God Help Me
    Map2Eigen* V,             // output - prim var
    const RowMajorMatrixXd& nix, // input - norm vec i dir x comp
    const RowMajorMatrixXd& niy, // input - norm vec i dir y comp
    const RowMajorMatrixXd& njx, // input - norm vec j dir x comp
    const RowMajorMatrixXd& njy, // input - norm vec j dir y comp
    RowMajorMatrixXd& T,        // input - grab the temperature 
    constants C               // input - constants C for case
    )
{
  int ni = V->Q[rhoid].cols()-2*C.num_ghost;
  int nj = V->Q[rhoid].rows()-2*C.num_ghost;

  int Lower_Begin[2]  = {0 , 0};
  int Lower_End[2]    = {0,ni-1};

  int Upper_Begin[2]  = {nj-1, 0};
  int Upper_End[2]    = { nj-1,ni-1};

  int Left_Begin[2]  = {0 , 0};
  int Left_End[2]    = {nj-1,0};

  int Right_Begin[2]  = {0, ni-1};
  int Right_End[2]    = { nj-1,ni-1};

  slipwallBC(V, Upper_Begin, Upper_End, nix, niy, njx, njy, T, C);
  slipwallBC(V, Lower_Begin, Lower_End, nix, niy, njx, njy, T, C);

  slipwallBC(V, Left_Begin, Left_End, nix, niy, njx, njy, T, C);
  slipwallBC(V, Right_Begin, Right_End, nix, niy, njx, njy, T, C);

}



void initialize(
    // This function initializes all of the primitive variable and it c
    // hecks for the correct case to initialize with. 
    Map2Eigen* V,           // output - Primitive variables 
    const RowMajorMatrixXd& xc_g,  // input - xcoord with ghosts
    const RowMajorMatrixXd& yc_g,  // input - ycoord with ghosts
    constants C            // input - constants for case number
    )
{

  int ni = xc_g.cols(); // n cells in i dir
  int nj = xc_g.rows();  // n cells in j dir
  // check for the same sizes
  if (V->Q[rhoid].cols() != ni || V->Q[rhoid].rows() != nj)
  {
    cerr << "ERROR: Dimension Mismatch in Initialization?!?!" << endl;
    exit(1);
  }
  // create important constants specified from papers
  double x, y, r, thet;
  double rhoL = 1.0;
  double rhoH = 2.0;

  double p0 = 2.5;
  double g = ACCEL;//0.5;

  double vPert = 0.0;
  double rmidd= 0.375;
  double circ = 2.0*PI*rmidd;
  double Bx = 0.0125;//125;
  double lambda = circ/20;
  double randPert;

  for (int j = 0; j < nj; j++)
  {
    for (int i = 0; i < ni; i++)
    {
      x = xc_g(j, i);
      y = yc_g(j, i);

      r = sqrt(x*x+y*y);
      thet = atan2(y,x);


      randPert = ((double) rand() / (RAND_MAX));
      if (r > rmidd*(1.0+0.02*randPert))//cos(1*abs(randPert)*thet/PI + randPert*PI)))
      {
        V->Q[rhoid](j,i) = rhoH;
      }
      else
      {
        V->Q[rhoid](j,i) = rhoL;
      }


      vPert =vPert;
      double decay =exp(-20.0*(r-rmidd)*(r-rmidd)/(0.25*rmidd*rmidd));
      double thetVar = randPert*(2.0*PI/lambda)*thet+randPert*PI;

      V->Q[uid](j,i) = vPert*decay;
      V->Q[vid](j,i) = vPert*decay;
      V->Q[wid](j,i) = 0.0; 
      V->Q[pid](j,i) = p0 - V->Q[rhoid](j,i)*g*(r-rmidd);
      V->Q[bxid](j,i) = Bx*sin(thet); 
      V->Q[byid](j,i) = -1*Bx*cos(thet); 
      V->Q[bzid](j,i) = 0.0; 
    }
  }
}




void computeVolume(
    // This function computes the volumes of the interior cell points. 
    RowMajorMatrixXd& Volume,      // output - Volume of the cells
    const RowMajorMatrixXd& xn,    // input - x coord for nodal vals
    const RowMajorMatrixXd& yn     // input - y coord for nodal vals
    )
{
  double dx1, dx2, dy1, dy2; // find the diagonals on the cell
  double cross; // cross product
  for (int i = 0; i < Volume.cols(); i++)
  {
    for (int j = 0; j < Volume.rows(); j++)
    {
      dx1 = xn(j,i) - xn(j+1,i+1); // bot L - top R
      dy1 = yn(j,i) - yn(j+1,i+1); // bot L - top R
      dx2 = xn(j+1,i) - xn(j,i+1); // top L - bot R
      dy2 = yn(j+1,i) - yn(j,i+1); // top L - bot R
      cross = abs(dx1*dy2 - dx2*dy1); // cross product
      Volume(j,i) = 0.5*cross; // assume w=1 volume of cell
    }
  } 
}

void computeNormalVectors(
    // This function takes in the previously calculated areas and
    // Finds the fraction in the x and y physical directions
    RowMajorMatrixXd& nix,    // output - i with A of normal vector in x phys
    RowMajorMatrixXd& niy,    // output - i comp but the yhat phys
    RowMajorMatrixXd& njx,    // output - j with A of normal vector in x phys
    RowMajorMatrixXd& njy,    // output - j dir but with yhat phys
    const RowMajorMatrixXd& xn,    // input - nodal x coordinate
    const RowMajorMatrixXd& yn,    // input - nodal y coordinate
    const RowMajorMatrixXd& Ai,    // input - Area pointed in i (x square grid)
    const RowMajorMatrixXd& Aj     // input - Area pointed in j (y square grid)
    )
{
  // NOTE i is the columns and j are the rows for A(row,col)
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      // See notes Sec 6 slide 75
      nix(j,i) = (yn(j+1,i) - yn(j,i)) / Ai(j,i); 
      niy(j,i) = -(xn(j+1,i) - xn(j,i)) / Ai(j,i);
    }
  }
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      njx(j,i) = -(yn(j,i+1) - yn(j,i)) / Aj(j,i);
      njy(j,i) = (xn(j,i+1) - xn(j,i)) / Aj(j,i);
    }
  }
}



void computeArea(
    // This function computes the area based on the node locations 
    // for both the i facing direction and the j facing direction
    RowMajorMatrixXd& Ai,          // output - i-direction areas (c, n)size    
    RowMajorMatrixXd& Aj,          // output - j-direction areas (n, c)size
    const RowMajorMatrixXd& xn,    // input - x coord for nodes 
    const RowMajorMatrixXd& yn     // input - y coord for nodes
    )
{
  // loop for Ai
  // NOTE i is the columns and j is the rows
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      double dx = xn(j,i) - xn(j+1,i); // compute bottom - top x
      double dy = yn(j,i) - yn(j+1,i); // compute bottom - top y
      Ai(j,i) = sqrt( dx*dx + dy*dy ); // pythag
    }
  }
  // loop for Aj
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      double dx = xn(j,i) - xn(j,i+1); // compute left - right x
      double dy = yn(j,i) - yn(j,i+1); // compute left - right y
      Aj(j,i) = sqrt( dx*dx + dy*dy);  // pytha
    }
  }
  // NOTE area is the dist times a w into the z dir which is just 1 for 2d
}











void extrapCopyCoords(
    // This function will copy the interior coord of the cells, but
    // will also extrapolate those coords to the ghost cells
    RowMajorMatrixXd& xc_g,        // output - cell x coord with ghosts 
    RowMajorMatrixXd& yc_g,        // output - cell y coord with ghosts
    const RowMajorMatrixXd& xc,    // input - cell x coord without ghosts
    const RowMajorMatrixXd& yc,    // input - cell y coord without ghosts
    constants C            // input - constants C for num ghost
    )
{
  int ni = xc.cols();
  int nj = xc.rows();
  // copy the interior values for both x and y
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      int i_g = i+C.num_ghost; // set iter to 3+i
      int j_g = j+C.num_ghost; // set iter to 3+j
      xc_g(j_g, i_g) = xc(j,i);
      yc_g(j_g, i_g) = yc(j,i);
    }
  }
  int bot = C.num_ghost - 1; // 1st ghost cell based off 
  int top = C.num_ghost + nj; // top index start
  int left = C.num_ghost -1; // same as bot
  int right = C.num_ghost + ni; // 1st ghost on right
  // extrapolate to the ghost cells now

  // j dir which means the interior need extrap
  for (int i = C.num_ghost; i < ni+C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bottom index - j is the filling of ghost
      xc_g(bot-j,i) = 2.0*xc_g(bot-j+1,i) - xc_g(bot-j+2,i);
      yc_g(bot-j,i) = 2.0*yc_g(bot-j+1,i) - yc_g(bot-j+2,i);
      // top index top + j is filling ghost
      xc_g(top+j,i) = 2.0*xc_g(top+j-1,i) - xc_g(top+j-2,i);
      yc_g(top+j,i) = 2.0*yc_g(top+j-1,i) - yc_g(top+j-2,i);
    }
  }
  // i dir which means the interior need extrap
  for (int i = 0; i < C.num_ghost; i ++)
  {
    for (int j = C.num_ghost; j < C.num_ghost+nj; j++)
    {
      // left index - i is the filling of ghost
      xc_g(j,left-i) = 2.0*xc_g(j,left-i+1) - xc_g(j,left-i+2);
      yc_g(j,left-i) = 2.0*yc_g(j,left-i+1) - yc_g(j,left-i+2);
      // right index + i is the filling of ghost
      xc_g(j,right+i) = 2.0*xc_g(j,right+i-1) - xc_g(j,right+i-2);
      yc_g(j,right+i) = 2.0*yc_g(j,right+i-1) - yc_g(j,right+i-2);
    }
  }
  //set to 1 out the corners for sake of clarity.
  for (int i = 0; i < C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bot left
      xc_g(j,i) = 1.0; yc_g(j,i) = 1.0;
      // bot right
      xc_g(j,right+i) = 1.0; yc_g(j,right+i) = 1.0;
      // top right
      xc_g(top+j,right+i) = 1.0; yc_g(top+j,right+i) = 1.0;
      // top left
      xc_g(top+j,i) = 1.0; yc_g(top+j,i) = 1.0;
    }
  }
}

void inputMesh(
    RowMajorMatrixXd& xn,         // Output - x coord nodes
    RowMajorMatrixXd& yn,         // Output - y coord nodes
    RowMajorMatrixXd& xc,         // Output - x coord centers
    RowMajorMatrixXd& yc,         // Output - y coord centers
    const string mesh)    // input - mesh location
{
  ifstream infile;
  infile.open(mesh.c_str());

  if (!infile){
    cerr << "ERROR: Unable to open the mesh file" << endl;
    exit(1);
  }

  double readval;
  string line;
  int vals[3];
  for(int i = 0; i < 3; i++ )
  {
    getline(infile, line);
    istringstream streamA(line);
    streamA >> readval;
    vals[i] = readval;
  }
  int ntot = vals[0];
  int nx = vals[1];
  int ny = vals[2];

  //double dat[ntot*2];
  VectorXd dat(ntot*2);
  int i = 0;
  while (infile.good())
    while (getline(infile, line))
    {
      istringstream streamA(line);
      while(streamA >> readval)
        dat(i) = readval;
      i++;
    } 

  xn.resize(ny, nx); //number of rows, number of columns
  yn.resize(ny, nx); //number of rows, number of columns
  xc.resize(ny-1, nx-1);
  yc.resize(ny-1, nx-1);

  int offset = ntot*2;
  for(int j = 0; j < ny; j++)
    for(int i = 0; i < nx; i++)
    {
      int xInd = j*nx+i;
      int yInd = j*nx+i+nx*ny;
      xn(j, i) = dat[xInd];
      yn(j, i) = dat[yInd];
    }
  for (int j = 0; j < xc.rows(); j++)
    for (int i = 0 ; i < xc.cols(); i++)
    {
      // take average of all four xvals from corn
      // order: botL + topL + botR + topR 
      xc(j,i) = 0.25*(xn(j,i) + xn(j+1,i) + 
          xn(j,i+1) + xn(j+1,i+1));
      yc(j,i) = 0.25*(yn(j,i) + yn(j+1,i) + 
          yn(j,i+1) + yn(j+1,i+1));
    }
}

void outputArrayMap(
    // This function will output any eigen matrix into a file
    string Address,        // input - folder location
    string FileName,       // input - File
    const Map2Eigen* out,         //input - matrix
    int n                  //input - iteration
    )
{
  string outString = FileName + "_rho";
  RowMajorMatrixXd outSingle = out->Q[rhoid];

  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[uid]; outString = FileName+"_u" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[vid]; outString = FileName+"_v" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[pid]; outString = FileName+"_p" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[bxid]; outString = FileName+"_bx" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[byid]; outString = FileName+"_by" ;
  outputArray(Address, outString, outSingle, n);
}

void outputArray(
    // This function will output any eigen matrix into a file
    string Address,        // input - folder location
    string FileName,       // input - File
    RowMajorMatrixXd& out,         //input - matrix
    int n                  //input - iteration
    )
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Suffix = ".txt";
  string Sub = "-";
  Address = Address + "/" + FileName + Sub + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component

  // Alot of options for outputting this is just nice to see 
  // where lines end and start
  //IOFormat CleanFmt(14,0,", ", "\n", "[", "]");
  IOFormat CleanFmt(14);
  outfile << out.format(CleanFmt) << endl; // output trans
  //outfile << setprecision(14) << out << endl; // output trans
}

void primToCons(
    Map2Eigen* U,           // output - Conserved vars
    const Map2Eigen* V)         // input - prim vars
{
  int ni = U->Q[rhoid].cols();
  int nj = U->Q[rhoid].rows();

  for(int j = 0; j < nj; j++)
    for(int i = 0; i < ni; i++)
    {
      U->Q[rhoid](j,i) = V->Q[rhoid](j,i); //U1 = V1
      U->Q[uid](j,i) = V->Q[rhoid](j,i)*(V->Q[uid](j,i));
      U->Q[vid](j,i) = V->Q[rhoid](j,i)*(V->Q[vid](j,i));
      U->Q[wid](j,i) = V->Q[rhoid](j,i)*(V->Q[wid](j,i));
      U->Q[pid](j,i) = V->Q[pid](j,i)/(GAMMA - 1.0) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[uid](j,i))*(V->Q[uid](j,i))) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[vid](j,i))*(V->Q[vid](j,i))) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[wid](j,i))*(V->Q[wid](j,i))) 
        + 0.5*V->Q[bxid](j,i)*((V->Q[bxid](j,i)))/MU0 
        + 0.5*V->Q[byid](j,i)*(V->Q[byid](j,i))/MU0 
        + 0.5*V->Q[bzid](j,i)*(V->Q[bzid](j,i))/MU0;

      U->Q[bxid](j,i) = V->Q[bxid](j,i); 
      U->Q[byid](j,i) = V->Q[byid](j,i); 
      U->Q[bzid](j,i) = V->Q[bzid](j,i); 
    }
}



void consToPrim(
    Map2Eigen* V,           // output - Conserved vars
    const Map2Eigen* U)        // input - prim vars
{
  int ni = V->Q[rhoid].cols();
  int nj = V->Q[rhoid].rows();
  for(int j = 0; j < nj; j++)
    for(int i = 0; i < ni; i++)
    {
      V->Q[rhoid](j,i) = U->Q[rhoid](j,i); //U1 = V1
      V->Q[uid](j,i)= U->Q[uid](j,i)/(U->Q[rhoid](j,i));
      V->Q[vid](j,i)= U->Q[vid](j,i)/(U->Q[rhoid](j,i));
      V->Q[wid](j,i)= U->Q[wid](j,i)/(U->Q[rhoid](j,i));
      V->Q[pid](j,i)= (GAMMA - 1.0) *(U->Q[pid](j,i)  
          - 0.5*U->Q[uid](j,i)*((U->Q[uid](j,i))/(U->Q[rhoid](j,i))) 
          - 0.5*U->Q[vid](j,i)*((U->Q[vid](j,i))/(U->Q[rhoid](j,i))) 
          - 0.5*U->Q[wid](j,i)*((U->Q[wid](j,i))/(U->Q[rhoid](j,i))) 
          - 0.5*U->Q[bxid](j,i)*(U->Q[bxid](j,i))/MU0 
          - 0.5*U->Q[byid](j,i)*(U->Q[byid](j,i))/MU0 
          - 0.5*U->Q[bzid](j,i)*(U->Q[bzid](j,i))/MU0);

      V->Q[bxid](j,i) = U->Q[bxid](j,i); 
      V->Q[byid](j,i) = U->Q[byid](j,i); 
      V->Q[bzid](j,i) = U->Q[bzid](j,i); 
    }
}

void cons2Prim(
    double V[],           // Output: Primitive variables
    double U[],           // Input: Conservative variables
    int size)             // Input: size of array must be evenly dividable by NEQ
{
  int index; 
  double U0[NEQ];
  for(int i = 0; i < size/NEQ; i++)
  {
    index = i*NEQ;
    for(int eq = 0; eq < NEQ; eq++)
      U0[eq] = U[index+eq], V[index+eq] = U0[eq];

    V[index+uid] = U0[uid]/U0[rhoid]; //rhou/rho
    V[index+vid] = U0[vid]/U0[rhoid]; //rhou/rho
    V[index+wid] = U0[wid]/U0[rhoid]; //rhou/rho
    V[index+pid] = (GAMMA-1.0)*(U0[pid] - 0.5*(U0[uid]*U0[uid] + U0[vid]*U0[vid] + U0[wid]*U0[wid])/U0[rhoid] - 0.5*((U0[bxid]*U0[bxid] + U0[byid]*U0[byid] + U0[bzid]*U0[bzid]))/MU0);
  }
}

//////////////////////////////// OLD //////////////////////


void computeFluxVL(double* FFLUX, double* ul, double* ur, double nxhat, double nyhat)
{
  // define constants
  double M_L, M_R, M_p, M_n;
  double beta_L, beta_R;
  double alpha_p, alpha_n;
  double c_p, c_n;
  double a_L, a_R;
  double U_L, U_R;
  double D_p, D_n;
  double p_bar_p, p_bar_n;
  double ht_L, ht_R; 
  double Fc[NEQ], Fp[NEQ];

  // grab the correct flux dimensions

  double vl[NEQ];
  double vr[NEQ];


  cons2Prim(vl, ul, NEQ);
  cons2Prim(vr, ur, NEQ);
  // compute sound speed
  a_L = sqrt(GAMMA*vl[pid] / vl[rhoid]); 
  a_R = sqrt(GAMMA*vr[pid] / vr[rhoid]); 

  // compute speed
  U_L = vl[uid]*nxhat + vl[vid]*nyhat;
  U_R = vr[uid]*nxhat + vr[vid]*nyhat;

  // compute mach # from speed Mach is a scalar***
  M_L = U_L/a_L;
  M_R = U_R/a_R;

  // Compute M_+ and M_- see notes sec06 slides > 74
  M_p = 0.25*(M_L + 1)*(M_L + 1);
  M_n = -0.25*(M_R - 1)*(M_R - 1);

  // Compute beta_LR for use in det flux cont
  beta_L = -mymax(0, 1 - int(abs(M_L)));
  beta_R = -mymax(0, 1 - int(abs(M_R)));

  // NOTE SIGN has specific inputs applies sign of b to a SIGN(a,b)
  alpha_p = 0.5*(1 + SIGN(1, M_L));
  alpha_n = 0.5*(1 - SIGN(1, M_R));

  c_p = alpha_p*(1+beta_L)*M_L - beta_L*M_p;
  c_n = alpha_n*(1+beta_R)*M_R - beta_R*M_n;

  // compute avg pressure
  p_bar_p = M_p*(-M_L + 2);
  p_bar_n = M_n*(-M_R - 2);

  // compute diff pos and neg
  D_p = alpha_p*(1+beta_L) - beta_L*p_bar_p;
  D_n = alpha_n*(1+beta_R) - beta_R*p_bar_n;

  // left and right enthalpy
  ht_L = (GAMMA/(GAMMA-1))*vl[pid]/vl[rhoid] 
    + 0.5*( vl[uid]*vl[uid] + 
        vl[vid]*vl[vid]);
  ht_R = (GAMMA/(GAMMA-1))*vr[pid]/vr[rhoid] 
    + 0.5*( vr[uid]*vr[uid] + 
        vr[vid]*vr[vid]);


  // compute convective flux contribution
  Fc[rhoid] = vl[rhoid]*a_L*c_p + vr[rhoid]*a_R*c_n;

  Fc[uid] = vl[rhoid]*a_L*c_p*vl[uid] 
    + vr[rhoid]*a_R*c_n*vr[uid];

  Fc[vid] = vl[rhoid]*a_L*c_p*vl[vid] 
    + vr[rhoid]*a_R*c_n*vr[vid];

  Fc[pid] = vl[rhoid]*a_L*c_p*ht_L 
    + vr[rhoid]*a_R*c_n*ht_R;

  // compute pressure flux contribution
  Fp[rhoid] = 0.0;
  Fp[uid] = D_p*nxhat*vl[pid] + 
    D_n*nxhat*vr[pid];
  Fp[vid] = D_p*nyhat*vl[pid] + 
    D_n*nyhat*vr[pid];
  Fp[pid] = 0.0;

  // loop over and add the two
  for (int eq = 0; eq < NEQ; eq++)
    FFLUX[eq] = Fc[eq] + Fp[eq];
}


void computeFluxRoe(double* FFLUX, double* ul, double* ur, double nxhat, double nyhat)
{
  // define roe averaged vars
  double rho_Roe;
  double u_Roe;
  double v_Roe;
  double U_hat_Roe_R, U_hat_Roe_L;
  double U_hat_Roe;
  double ht_L, ht_R;
  double ht_Roe;
  double a_Roe;
  double R_Roe;

  // Get prim var
  double vl[NEQ];
  double vr[NEQ];

  cons2Prim(vl, ul, NEQ);
  cons2Prim(vr, ur, NEQ);

  // need delta's to compute the characteristic vars
  double drho, du, dv, dp;
  double dw1, dw2, dw3, dw4;

  // eigen vals
  double lambda1, lambda2, lambda3, lambda4;
  double eps = 0.1; // correction term
  // right eigen vectors
  double r1[NEQ], r2[NEQ], r3[NEQ], r4[NEQ];

  // Used to find the 1st order contribution to the flux determination
  double F_L[NEQ];
  double F_R[NEQ];

  // 2nd order contribution
  double sum2ndOrder;

  //int ni = vl[uid];
  //int nj = vl[uid].rows();

  // compute roe avgd quantities like before
  R_Roe = sqrt(vr[rhoid]/vl[rhoid]);
  rho_Roe = R_Roe*vl[rhoid];
  u_Roe = (R_Roe*vr[uid] + vl[uid]) / (R_Roe + 1);
  v_Roe = (R_Roe*vr[vid] + vl[vid]) / (R_Roe + 1);
  // Find the 2D effect on the speed
  U_hat_Roe = u_Roe*nxhat + v_Roe*nyhat;// get the speed

  // Roe averaged vars
  ht_L = (GAMMA/(GAMMA-1))*vl[pid]/vl[rhoid] 
    + 0.5*( vl[uid]*vl[uid] + vl[vid]*vl[vid]);
  ht_R = (GAMMA/(GAMMA-1))*vr[pid]/vr[rhoid] 
    + 0.5*( vr[uid]*vr[uid] + vr[vid]*vr[vid]);
  ht_Roe = (R_Roe*ht_R + ht_L) / (R_Roe + 1);

  // sound speed with the energy uu+vv (2D)
  a_Roe = sqrt((GAMMA-1)*( ht_Roe - 0.5*(u_Roe*u_Roe + v_Roe*v_Roe)));

  // Two of the eigen values are repeated
  lambda1 = abs(U_hat_Roe); // speed 
  lambda2 = abs(U_hat_Roe); // speed
  lambda3 = abs(U_hat_Roe + a_Roe); // u+a 1D
  lambda4 = abs(U_hat_Roe - a_Roe); // u-a 1D

  // apply expansion fan fix
  // fix the repeated eigen values (u) 
  if (lambda1 < 2*eps*a_Roe)
  {
    lambda1 = lambda1*lambda1/(4*eps*a_Roe) + eps*a_Roe;
    lambda2 = lambda2*lambda2/(4*eps*a_Roe) + eps*a_Roe;
  }
  // fix the u+a and u-a
  if (lambda3 < 2*eps*a_Roe)
    lambda3 = lambda3*lambda3/(4*eps*a_Roe) + eps*a_Roe;
  if (lambda4 < 2*eps*a_Roe)
    lambda4 = lambda4*lambda4/(4*eps*a_Roe) + eps*a_Roe;

  // fill in the right eigen vectors
  r1[rhoid] = 1.0;
  r1[uid] = u_Roe;
  r1[vid] = v_Roe;
  r1[pid] = 0.5*(u_Roe*u_Roe + v_Roe*v_Roe);

  r2[rhoid] = 0.0;
  r2[uid] = rho_Roe * nyhat;
  r2[vid] = -rho_Roe * nxhat;
  r2[pid] = rho_Roe*(u_Roe*nyhat - v_Roe*nxhat);

  r3[rhoid] = (0.5*rho_Roe/a_Roe);
  r3[uid] = (0.5*rho_Roe/a_Roe) * ( u_Roe + a_Roe*nxhat );
  r3[vid] = (0.5*rho_Roe/a_Roe) * ( v_Roe + a_Roe*nyhat );
  r3[pid] = (0.5*rho_Roe/a_Roe) * ( ht_Roe + a_Roe*U_hat_Roe );

  r4[rhoid] = (-0.5*rho_Roe/a_Roe);
  r4[uid] = (-0.5*rho_Roe/a_Roe) * ( u_Roe - a_Roe*nxhat );
  r4[vid] = (-0.5*rho_Roe/a_Roe) * ( v_Roe - a_Roe*nyhat );
  r4[pid] = (-0.5*rho_Roe/a_Roe) * ( ht_Roe - a_Roe*U_hat_Roe );

  // compute delta's
  drho = vr[rhoid] - vl[rhoid];
  du = vr[uid] - vl[uid];
  dv = vr[vid] - vl[vid];
  dp = vr[pid] - vl[pid];

  // compute characteristic variables
  dw1 = drho - dp/(a_Roe*a_Roe);
  dw2 = du*nyhat - dv*nxhat;
  dw3 = du*nxhat + dv*nyhat + dp/(rho_Roe*a_Roe);
  dw4 = du*nxhat + dv*nyhat - dp/(rho_Roe*a_Roe);

  // convert the reimann problem to 1d
  U_hat_Roe_L = ( vl[uid]*nxhat + vl[vid]*nyhat);
  U_hat_Roe_R = ( vr[uid]*nxhat + vr[vid]*nyhat);

  // grab the first order contributions
  F_L[rhoid] = vl[rhoid]*U_hat_Roe_L;
  F_R[rhoid] = vr[rhoid]*U_hat_Roe_R;

  F_L[uid] = vl[rhoid]*vl[uid]*U_hat_Roe_L + vl[pid]*nxhat;
  F_R[uid] = vr[rhoid]*vr[uid]*U_hat_Roe_R + vr[pid]*nxhat;

  F_L[vid] = vl[rhoid]*vl[vid]*U_hat_Roe_L + vl[pid]*nyhat;
  F_R[vid] = vr[rhoid]*vr[vid]*U_hat_Roe_R + vr[pid]*nyhat;

  F_L[pid] = vl[rhoid]*ht_L*U_hat_Roe_L;
  F_R[pid] = vr[rhoid]*ht_R*U_hat_Roe_R;

  // Combine 1st order and 2nd order
  // sum over the right vectors
  for (int eq = 0; eq < NEQ; eq++)
  {
    sum2ndOrder=0.5*(abs(lambda1)*dw1*r1[eq]
        + abs(lambda2)*dw2*r2[eq]
        + abs(lambda3)*dw3*r3[eq]
        + abs(lambda4)*dw4*r4[eq]);
    // flux = 1st order - 2ndordersum(vec(r)*vec(dw)*abs(lambda)
    FFLUX[eq] = 0.5*(F_L[eq]+F_R[eq]) - sum2ndOrder;
  }
}


void computeFluxHLLD(double F[], double UL[], double UR[], double& nxhat, double& nyhat, int ForG)
{


  double CL = computeMaxSpeed(UL, ForG);
  double CR = computeMaxSpeed(UR, ForG);

  double rhoL = UL[rhoid];
  double rhoR = UR[rhoid];
  double pL = UL[pid];
  double pR = UR[pid];



  double uL = (UL[uid]*nxhat + UL[vid]*nyhat)/rhoL; 
  double uR = (UR[uid]*nxhat + UR[vid]*nyhat)/rhoR; 

  // avg the mag
  double bparrL = (UL[bxid]*nxhat + UL[byid]*nyhat); 
  double bparrR = (UR[bxid]*nxhat + UR[byid]*nyhat); 
  double bparr = (bparrL+bparrR)/2; 

  double byL = nxhat*(UL[byid]*nxhat - UL[bxid]*nyhat); 
  double byR = nxhat*(UR[byid]*nxhat - UR[bxid]*nyhat); 

  double bzL = UL[bzid]*(nxhat*nxhat+nyhat*nyhat);
  double bzR = UR[bzid]*(nxhat*nxhat+nyhat*nyhat);

  double vL = nxhat*(UL[vid]*nxhat - UL[uid]*nyhat); 
  double vR = nxhat*(UR[vid]*nxhat - UR[uid]*nyhat); 

  double wL = UL[wid]*(nxhat*nxhat+nyhat*nyhat);
  double wR = UR[wid]*(nxhat*nxhat+nyhat*nyhat);

  double SL = mymin(uL, uR) - mymax(CL, CR);
  double SR = mymax(uL, uR) + mymax(CL, CR);

  // calculate the pressure across the jump then use this for the SM
  double pTL = pL + 0.5*(UL[bxid]*UL[bxid]+UL[byid]*UL[byid]+UL[bzid]*UL[bzid])/MU0;
  double pTR = pR + 0.5*(UR[bxid]*UR[bxid]+UR[byid]*UR[byid]+UR[bzid]*UR[bzid])/MU0;
  double pTs = ((SR-uR)*rhoR*pTL - (SL-uL)*rhoL*pTR + rhoL*rhoR*(SR-uR)*(SL-uL)*(uR-uL))/
    ((SR-uR)*rhoR-(SL-uL)*rhoL);

  // can get uL and uR etc
  double SM = ((SR-uR)*rhoR*uR-(SL-uL)*rhoL*uL-pTR+pTL)/((SR-uR)*rhoR-(SL-uL)*rhoL);


  // compute intermediate states
  double uLs = SM;
  double uRs = SM;

  double pTLs = pTs;
  double pTRs = pTs;

  double rhoLs = rhoL*(SL-uL)/(SL-SM);
  double rhoRs = rhoR*(SR-uR)/(SR-SM);

  double vLs = vL-bparr*byL*(SM-uL)/(rhoL*(SL-uL)*(SL-SM)-bparr*bparr);
  double vRs = vR-bparr*byR*(SM-uR)/(rhoR*(SR-uR)*(SR-SM)-bparr*bparr);

  double byLs = byL*(rhoL*(SL-uL)*(SL-uL) - bparr*bparr)/(rhoL*(SL-uL)*(SL-SM)-bparr*bparr);
  double byRs = byR*(rhoR*(SR-uR)*(SR-uR) - bparr*bparr)/(rhoR*(SR-uR)*(SR-SM)-bparr*bparr);

  double wLs = wL - bparr*bzL*((SM-uL)/(rhoL*(SL-uL)*(SL-SM)-bparr*bparr));
  double wRs = wR - bparr*bzR*((SM-uR)/(rhoR*(SR-uR)*(SR-SM)-bparr*bparr));

  double bzLs = bzL*((rhoL*(SL-uL)*(SL-uL)-bparr*bparr)/(rhoL*(SL-uL)*(SL-SM)-bparr*bparr));
  double bzRs = bzR*((rhoR*(SR-uR)*(SR-uR)-bparr*bparr)/(rhoR*(SR-uR)*(SR-SM)-bparr*bparr));

  double eLs = ((SL-uL)*UL[pid] - pTL*uL + pTs*SM + bparr*((vL*byL+wL*bzL) - (uLs*bparr+vLs*byLs+wLs*bzLs)))/(SL-SM);
  double eRs = ((SR-uR)*UR[pid] - pTR*uR + pTs*SM + bparr*((vR*byR+wR*bzR) - (uRs*bparr+vRs*byRs+wRs*bzRs)))/(SR-SM);

  // compute next state

  double uLss = SM;
  double uRss = SM;
  double pTRss = pTs;
  double pTLss = pTs;

  double rhoLss = rhoLs;
  double rhoRss = rhoRs;

  double SLs = SM-abs(bparr)/sqrt(rhoLs);
  double SRs = SM+abs(bparr)/sqrt(rhoRs);

  double vLss = (sqrt(rhoLs)*vLs+sqrt(rhoRs)*vRs+(byRs-byLs)*SIGN(1,bparr))/(sqrt(rhoLs)+sqrt(rhoRs));
  double vRss = vLss;

  double wLss = (sqrt(rhoLs)*wLs+sqrt(rhoRs)*wRs+(bzRs-bzLs)*SIGN(1,bparr))/(sqrt(rhoLs)+sqrt(rhoRs));
  double wRss = wLss;

  double byLss = (sqrt(rhoLs)*byRs+sqrt(rhoRs)*byLs+sqrt(rhoLs*rhoRs)*(vRs-vLs)*SIGN(1,bparr))/(sqrt(rhoLs)+sqrt(rhoRs));
  double byRss = byLss;


  double bzLss = (sqrt(rhoLs)*bzRs+sqrt(rhoRs)*bzLs+sqrt(rhoLs*rhoRs)*(wRs-wLs)*SIGN(1,bparr))/(sqrt(rhoLs)+sqrt(rhoRs));
  double bzRss = bzLss;

  double eLss = eLs - sqrt(rhoLs)*((uLs*bparr+vLs*byLs+wLs*bzLs) - (uLss*bparr+vLss*byLss+wLss*bzLss));
  double eRss = eRs - sqrt(rhoRs)*((uRs*bparr+vRs*byRs+wRs*bzRs) - (uRss*bparr+vRss*byRss+wRss*bzRss));

  double ULs[NEQ], URs[NEQ];
  double ULss[NEQ], URss[NEQ];

  ULs[rhoid] = rhoLs;
  ULs[uid] = uLs;
  ULs[vid] = vLs;
  ULs[wid] = wLs;
  ULs[pid] = eLs;
  ULs[bxid] = bparr;
  ULs[byid] = byLs;
  ULs[bzid] = bzLs;

  ULss[rhoid] = rhoLss;
  ULss[uid] =     uLss;
  ULss[vid] =     vLss;
  ULss[wid] =     wLss;
  ULss[pid] =     eLss;
  ULss[bxid] =  bparr;
  ULss[byid] =   byLss;
  ULss[bzid] =   bzLss;

  URs[rhoid] = rhoRs;
  URs[uid] = uRs;
  URs[vid] = vRs;
  URs[wid] = wRs;
  URs[pid] = eRs;
  URs[bxid] = bparr;
  URs[byid] = byRs;
  URs[bzid] = bzRs;

  URss[rhoid] = rhoRss;
  URss[uid] =     uRss;
  URss[vid] =     vRss;
  URss[wid] =     wRss;
  URss[pid] =     eRss;
  URss[bxid] =  bparr;
  URss[byid] =   byRss;
  URss[bzid] =   bzRss;

  if (ForG == 0)
  {
    if (SL > 0)
      fFlux(F, UL, nxhat, nyhat);
    if (SL <= 0 && SLs >= 0)
      fFlux(F, ULs, nxhat, nyhat);
    if (SLs <= 0 && SM >= 0)
      fFlux(F, ULss, nxhat, nyhat);
    if (SM <= 0 && SRs >= 0)
      fFlux(F, URss, nxhat, nyhat);
    if (SRs <= 0 && SR >= 0)
      fFlux(F, URs, nxhat, nyhat);
    if (SR < 0)
      fFlux(F, UR, nxhat, nyhat);
  }
  else
  {
    if (SL > 0)
      gFlux(F, UL, nxhat, nyhat);
    if (SL <= 0 && SLs >= 0)
      gFlux(F, ULs, nxhat, nyhat);
    if (SLs <= 0 && SM >= 0)
      gFlux(F, ULss, nxhat, nyhat);
    if (SM <= 0 && SRs >= 0)
      gFlux(F, URss, nxhat, nyhat);
    if (SRs <= 0 && SR >= 0)
      gFlux(F, URs, nxhat, nyhat);
    if (SR < 0)
      gFlux(F, UR, nxhat, nyhat);
  }

  /*
     if (ForG == 0)
     {
     fFlux(FA, UA);
     fFlux(FB, UB);
     }
     else if (ForG == 1)
     {
     gFlux(FA, UA);
     gFlux(FB, UB);
     }
     else
     {
     cerr << "ERROR: Wrong index for direction" << endl;
     exit(-1);
     }

  //for(int eq=0; eq < NEQ; eq++)
  if (lambdaA > 0)
  for(int eq = 0; eq < NEQ; eq++)
  F[eq] = FA[eq];
  else if (lambdaB < 0)
  for(int eq = 0; eq < NEQ; eq++)
  F[eq] = FB[eq];
  else if (lambdaB >=0 && lambdaA <=0)
  for(int eq = 0; eq < NEQ; eq++)
  F[eq] = (lambdaB*FA[eq] - lambdaA*FB[eq]+lambdaA*lambdaB*(UB[eq]-UA[eq]))/(lambdaB-lambdaA);
  else
  {
  cerr << "ERROR: HLL Flux Calculation Lambda Unphysical" << endl;
  exit(-1);
  }
  */
}

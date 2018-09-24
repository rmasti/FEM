#include "mhdRT.hpp"

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);

  constants C;
  C.f_limiter = 7;
  C.num_ghost = 3;
  C.cfl = 0.45;
  C.nmax = 10000;
  C.wint = 500;
  C.pint = 1;
  double A = (2-1)/(1.0+2.0);
  double tend = 6.0/sqrt(A*ACCEL*2);
  double dt;
  MatrixXd time(1, C.nmax);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  string mesh = "../mesh/debugMatlab.msh";
  string outputFolder = "./output/";


  RowMajorMatrixXd xcL_g, ycL_g;
  RowMajorMatrixXd nixL, niyL;
  RowMajorMatrixXd njxL, njyL;
  RowMajorMatrixXd AjL, AiL;
  RowMajorMatrixXd VolumeL;

  //meshBlock(mesh, outputFolder, C);
  // grab the information needed to run initialize and the fvm locally
  if(rank == 0)
    cout << " Grabbing Geometry " << endl;
  MPI_Comm com2d;
  com2d = meshBlock(mesh, outputFolder, xcL_g, ycL_g, nixL, niyL, njxL, njyL, AiL, AjL, VolumeL, C);

  int coordMax[2]={0, 0};
  MPI_Cart_coords(com2d, rank, 2, coordMax);

  MPI_Allreduce(&coordMax[0], &coordMax[0], 1, MPI_INT, MPI_MAX, com2d);
  MPI_Allreduce(&coordMax[1], &coordMax[1], 1, MPI_INT, MPI_MAX, com2d);

  MPI_Comm_rank(com2d, &rank);

  int nic_g = xcL_g.cols();
  int njc_g = ycL_g.rows();

  // Create the map2Eigen vars for storage with ghost layers
  Map2Eigen *V = new Map2Eigen(njc_g, nic_g, NEQ);
  Map2Eigen *U = new Map2Eigen(njc_g, nic_g, NEQ);
  Map2Eigen *U_RK = new Map2Eigen(njc_g, nic_g, NEQ);

  int njc = VolumeL.rows();
  int nic = VolumeL.cols();
  Map2Eigen *S    = new Map2Eigen(njc, nic, NEQ);
  Map2Eigen *Res  = new Map2Eigen(njc, nic, NEQ);
  //Direction dependent
  Map2Eigen *F    = new Map2Eigen(njc, nic+1, NEQ);
  Map2Eigen *U_L  = new Map2Eigen(njc, nic+1, NEQ);
  Map2Eigen *U_R  = new Map2Eigen(njc, nic+1, NEQ);

  Map2Eigen *G    = new Map2Eigen(njc+1, nic, NEQ);
  Map2Eigen *U_T  = new Map2Eigen(njc+1, nic, NEQ);
  Map2Eigen *U_B  = new Map2Eigen(njc+1, nic, NEQ);


  RowMajorMatrixXd xcL = xcL_g.block(C.num_ghost, C.num_ghost, njc, nic);
  RowMajorMatrixXd ycL = ycL_g.block(C.num_ghost, C.num_ghost, njc, nic);
  // initialize prim var

  if(rank == 0)
    cout << " Initializing " << endl;
 
  initialize(V, xcL_g, ycL_g, C);
  primToCons(U,V);

  U_RK = U;

  // send and recv func
  //setBcSend(U, nixL, niyL, njxL, njyL, com2d,  C);
  //setBcRecv(U, com2d, C);
  //MPI_Barrier(com2d);

  // recv while you send
  //cout << " setting bc " << endl;
  mpiSetBc(U, nixL, niyL, njxL, njyL, com2d,  C);

  MPI_Barrier(com2d);

  outputArrayMap(outputFolder, "UL", U, rank);
  outputArrayMap(outputFolder, "VL", V, rank);
  outputArray(outputFolder, "xcL_g", xcL_g, rank+size);

  dt = computeTimeStepU(VolumeL, AiL, AjL, nixL, niyL, njxL, njyL, U, C);
  MPI_Allreduce(&dt, &dt, 1, MPI_DOUBLE, MPI_MIN, com2d);

  int n=0;

  //cout << "Before stitch" << endl;

  MPI_Barrier(com2d);
  stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);

  if(rank == 0)
    cout << " Entering Time Loop " << endl;

  
  /*
  while(time(0, n) < tend)
  {
    for (int k = 0; k < RKORDER; k++)
    {
      setBcSend(U, nixL, niyL, njxL, njyL, com2d,  C);

      setBcRecv(U, com2d, C);
      MPI_Barrier(com2d);

      computeSourceTerm(S, U_RK, xcL, ycL, C);

      MUSCL(U_L, U_R, U_B, U_T, U_RK, C);

      compute2dFlux(F, G, U_L, U_R, U_B, U_T, njxL, njyL, nixL, niyL, C);

      computeRes(Res, S, F, G, AjL, AiL, VolumeL, C);

      rungeKutta(U_RK, U, Res, VolumeL, k, dt, C);
    }

    U=U_RK;


    //cout<< " Made Through RK LOOP" << endl;
    if(time(0,n)+dt > tend)// end at exact time
      dt = tend-time(0,n);
    else
    {
      dt = computeTimeStepU(VolumeL, AiL, AjL, nixL, niyL, njxL, njyL, U, C);
      MPI_Allreduce(&dt, &dt, 1, MPI_DOUBLE, MPI_MIN, com2d);
    }

    time(0,n+1) = time(0,n) + dt;
    n++;


    MPI_Comm_rank(com2d, &rank);
    if (rank == 0)
      if(n%C.pint == 0)
        cout << "time  = " << time(0,n) << ", n = " << n << endl;

    if(n%C.wint == 0) // Write to file by rank 0 gather and output
      stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);

  }

  stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);
  */
  MPI_Finalize();
  return 0;
}

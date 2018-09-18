#include "mhdRT.hpp"

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);

  constants C;
  C.f_limiter = 7;
  C.num_ghost = 3;
  C.cfl = 0.45;
  C.nmax = 10000;
  C.wint = 100;
  C.pint = 10;
  double A = (2-1)/(1.0+2.0);
  double tend = 6.0/sqrt(A*ACCEL*2);

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
  MPI_Comm com2d;
  com2d = meshBlock(mesh, outputFolder, xcL_g, ycL_g, nixL, niyL, njxL, njyL, AiL, AjL, VolumeL, C);
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

  // initialize prim var
  initialize(V, xcL_g, ycL_g, C);
  primToCons(U,V);

  U_RK = U;

  outputArrayMap(outputFolder, "U", U, rank);
  outputArrayMap(outputFolder, "V", V, rank);
  outputArray(outputFolder, "xcL_g", xcL_g, rank+size);

  setBcSend(U, nixL, niyL, njxL, njyL, com2d,  C);

/* 
  int coord[2];
  int u, l, d, r;
  MPI_Cart_coords(com2d, rank, 2, coord);
  MPI_Cart_shift(com2d, 0, 1, &d, &u);
  MPI_Cart_shift(com2d, 1, 1, &l, &r);
      
  printf("My rank is %d: My coordinates are %d, %d\nMy neighbors are left:%d right:%d up:%d down:%d\n", rank, coord[0], coord[1], l, r, u ,d);

  if(rank == 2)
  {
    cout << endl;
    cout << xcL_g << endl;
    cout << endl;
    cout << ycL_g << endl;
  }

  */

  MPI_Finalize();
  return 0;
}

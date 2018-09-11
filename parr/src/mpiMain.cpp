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

  meshBlock(mesh, outputFolder, C);



  MPI_Finalize();
  return 0;
}

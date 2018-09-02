# include "mhdRT.hpp"

int main (int argc, char * argv[])
{
  MPI_Init(& argc, &argv);


  int rank; int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size); 

  int ni = 3, nj = 5;


  if(rank == 0)
  {
    Map2Eigen *Var = new Map2Eigen(ni, nj, NEQ);

    for (int i = 0; i < ni*nj*NEQ; i++)
      Var->Q_raw[i] = i;

    Var->Q[rhoid](0,0) = 100;


    int dest = 1;
    int tag = 666;
    cout << "Sending from rank 0 process to rank 1 this array:" << endl;
    cout << Var->Q[0] << endl;
    double a[3] = {1.2, 3.4, 5.6};
    MPI_Send(Var->Q_raw, ni*nj*NEQ, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
    //MPI_Send(a, 3, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
  }

  if(rank == 1)
  {
    int tag = 666;

    Map2Eigen *VarIn = new Map2Eigen(ni, nj, NEQ);

    int source = 0;
    MPI_Status status;

    MPI_Recv(VarIn->Q_raw, ni*nj*NEQ, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
   
    //double aIn[3];
    //MPI_Recv(aIn, 3, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);

    cout << "Received this array from rank 0:" << endl;
    cout << VarIn->Q[0] << endl;

  }

  MPI_Finalize();

  return 0;

}

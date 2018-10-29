#include "mhdRT.hpp"


// constructor


Map2Eigen :: Map2Eigen(int nj, int ni, int nequ):
  Q_raw(new double[nj*ni*nequ]),
  Q{{Q_raw,   nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+1, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+2, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+3, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+4, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+5, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+6, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
    {Q_raw+7, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) } }
{

}




void stitchMap2EigenWrite(string Address, string FileName,const Map2Eigen* IN, const int n, const int *coord, MPI_Comm &com2d, constants C)
{
  int njc_g = IN->Q[rhoid].rows(); // get interiors
  int nic_g = IN->Q[rhoid].cols(); // get interiors


  int njc = njc_g-2*C.num_ghost; // get interiors
  int nic = nic_g-2*C.num_ghost; // get interiors

  int rank;
  int size;
  MPI_Comm_rank(com2d, &rank);
  MPI_Comm_size(com2d, &size);
  MPI_Status status;
  MPI_Request requestOut;
  MPI_Request requestIn;

  //printf("max j = %d; max i = %d\n", coord[0], coord[1]);

  int nip = coord[1]+1;
  int njp = coord[0]+1;

  int r;

  Map2Eigen *temp = new Map2Eigen(njc, nic, NEQ);
  int sizeDat = njc*nic*NEQ; 
  if (rank != 0) // send the data
  {
    for(int eq = 0; eq < NEQ; eq++)
      temp->Q[eq] = IN->Q[eq].block(C.num_ghost, C.num_ghost, njc, nic);
    MPI_Send(temp->Q_raw, sizeDat, MPI_DOUBLE, 0, rank*111, com2d);
  }

  //MPI_Wait(&requestOut,&status);
  //MPI_Barrier(com2d);
  if (rank == 0) // need to recv data
  {
    Map2Eigen *out = new Map2Eigen(njc*njp, nic*nip, NEQ);
    for(int j = 0; j < njp; j++)
      for(int i = 0; i < nip; i++)
        if(i == 0 && j == 0)// your bottom left fill as such
        {
          for(int eq = 0; eq < NEQ; eq++)
            out->Q[eq].block(0,0, njc, nic) = IN->Q[eq].block(C.num_ghost, C.num_ghost, njc, nic);
        }
        else
        {
          r = (j*nip)+i;

          MPI_Recv(temp->Q_raw, sizeDat, MPI_DOUBLE, r, 111*r, com2d, &status); 
          for(int eq = 0; eq < NEQ; eq++)
            out->Q[eq].block(njc*j,nic*i, njc, nic) = temp->Q[eq];
        }      
    outputArrayMap(Address, FileName, out, n);
  }
  delete temp->Q_raw; temp->Q_raw = NULL;

  MPI_Barrier(com2d);
}





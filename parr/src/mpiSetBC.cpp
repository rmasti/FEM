#include "mhdRT.hpp"

void setBcSend(Map2Eigen* U, const RowMajorMatrixXd& nix, const RowMajorMatrixXd& niy, const RowMajorMatrixXd& njx, const RowMajorMatrixXd& njy, MPI_Comm& com2d,  constants C)
{
  int rank;
  int size;
  MPI_Comm_rank(com2d, &rank);
  MPI_Comm_size(com2d, &size);

  int u, l, d, r;
  computeCoord(l, r, d, u, com2d);

  int njc = U->Q[rhoid].rows();
  int nic = U->Q[rhoid].cols();
  Map2Eigen *tempLR  = new Map2Eigen(njc , C.num_ghost, NEQ);

  if(rank < r && (l > r || l == r))// physical boundary on the left ergo I need to send my cells to the right
  {
  
    int icl = C.num_ghost;
    for(int i = 0; i < C.num_ghost; i++)
    {
      // no flip sig
      tempLR->Q[rhoid].col(i) = U->Q[rhoid].col(icl+i);
      tempLR->Q[wid].col(i) = U->Q[wid].col(icl+i);
      tempLR->Q[bzid].col(i) = U->Q[bzid].col(icl+i);
      tempLR->Q[pid].col(i) = U->Q[pid].col(icl+i);

      // flip sign
      tempLR->Q[uid].col(i) = -1.0*U->Q[uid].col(icl+i);
      tempLR->Q[vid].col(i) = -1.0*U->Q[vid].col(icl+i);
      tempLR->Q[bxid].col(i) = -1.0*U->Q[bxid].col(icl+i);
      tempLR->Q[byid].col(i) = -1.0*U->Q[byid].col(icl+i);
    }
    //MPI_Send(tempLR->Q_raw, njc*C.num_ghost, MPI_DOUBLE, r, 222, com2d);
  }


  printf("My rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
}



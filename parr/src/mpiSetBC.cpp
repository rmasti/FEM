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

  MPI_Request request; 
  MPI_Status status; 

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// SEND LEFT AND RIGHT DATA  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
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
    // take my data and send it to the processor on the left which is periodic phys boundary
    //MPI_Isend(tempLR->Q_raw, njc*C.num_ghost, MPI_DOUBLE, l, 222, com2d, &request); 
    MPI_Send(tempLR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, l, 222, com2d);
    //printf("My left on phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }
  else // not physical boundary send my stuff on left of my interior to the rank to the left
  {
    int icl = C.num_ghost;
    for(int i = 0; i < C.num_ghost; i++)
    {
      // no flip sign just grab my cell
      tempLR->Q[rhoid].col(i) = U->Q[rhoid].col(icl+i);
      tempLR->Q[wid].col(i) = U->Q[wid].col(icl+i);
      tempLR->Q[bzid].col(i) = U->Q[bzid].col(icl+i);
      tempLR->Q[pid].col(i) = U->Q[pid].col(icl+i);
      tempLR->Q[uid].col(i) = U->Q[uid].col(icl+i);
      tempLR->Q[vid].col(i) = U->Q[vid].col(icl+i);
      tempLR->Q[bxid].col(i) = U->Q[bxid].col(icl+i);
      tempLR->Q[byid].col(i) = U->Q[byid].col(icl+i);
    }
    // take my data and send it to the processor on the left which is not phys boundary
    //MPI_Isend(tempLR->Q_raw, njc*C.num_ghost, MPI_DOUBLE, l, 222, com2d, &request); // tag 222 is right
    MPI_Send(tempLR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, l, 222, com2d);
    //printf("My left not on phys my rank is %d: My neighbors are leftInt:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }

  // right physical boundary
  if(rank > r && (r < l || l == r))
  {

    int icr = (nic-C.num_ghost)-1; // first cell layer on right 
    for(int i = 0; i < C.num_ghost; i++)
    {
      // no flip sig
      tempLR->Q[rhoid].col(i) = U->Q[rhoid].col(icr-i);
      tempLR->Q[wid].col(i) = U->Q[wid].col(icr-i);
      tempLR->Q[bzid].col(i) = U->Q[bzid].col(icr-i);
      tempLR->Q[pid].col(i) = U->Q[pid].col(icr-i);

      // flip sign
      tempLR->Q[uid].col(i) = -1.0*U->Q[uid].col(icr-i);
      tempLR->Q[vid].col(i) = -1.0*U->Q[vid].col(icr-i);
      tempLR->Q[bxid].col(i) = -1.0*U->Q[bxid].col(icr-i);
      tempLR->Q[byid].col(i) = -1.0*U->Q[byid].col(icr-i);
    }
    // take my data and send it to the processor to the right of my which is periodic phys boundary
    //MPI_Isend(tempLR->Q_raw, njc*C.num_ghost, MPI_DOUBLE, r, 111, com2d, &request);// tag 111 is left
    MPI_Send(tempLR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, r, 111, com2d);
    //printf("My right on phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }
  else
  {
    int icr = (nic-C.num_ghost)-1; // first cell layer on right 
    for(int i = 0; i < C.num_ghost; i++)
    {
      // no flip sig
      tempLR->Q[rhoid].col(i) = U->Q[rhoid].col(icr-i);
      tempLR->Q[wid].col(i) = U->Q[wid].col(icr-i);
      tempLR->Q[bzid].col(i) = U->Q[bzid].col(icr-i);
      tempLR->Q[pid].col(i) = U->Q[pid].col(icr-i);
      tempLR->Q[uid].col(i) =  U->Q[uid].col(icr-i);
      tempLR->Q[vid].col(i) =  U->Q[vid].col(icr-i);
      tempLR->Q[bxid].col(i) = U->Q[bxid].col(icr-i);
      tempLR->Q[byid].col(i) = U->Q[byid].col(icr-i);
    }
    // take my data and send it to the processor to the right which is not phys boundary
    //MPI_Isend(tempLR->Q_raw, njc*C.num_ghost, MPI_DOUBLE, r, 111, com2d, &request);// tag 111 is left
    MPI_Send(tempLR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, r, 111, com2d);
    //printf("My right not on phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }

  //tempLR->Qraw -> Nu
  delete[] tempLR->Q_raw; tempLR->Q_raw = NULL;

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// SEND BOTT AND TOP  DATA  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  Map2Eigen *tempBT  = new Map2Eigen(C.num_ghost, nic, NEQ);

  if(d < 0)// physical boundary on bottom 
  {
    int jfb = 0;
    int jgb;
    int iff;

    int jcb;
    double uvel, vvel, bx, by, nx,ny;

    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {

        iff = i - C.num_ghost;// index for normvecs
        jcb = C.num_ghost+j; // cell interior
        jgb = (C.num_ghost-1) - j; // ghost cell row-layer

        // bottom boundary no penetration conducting wall
        nx = njx(jfb,iff);
        ny = njy(jfb,iff);
        uvel = U->Q[uid](jcb+j,i);
        vvel = U->Q[vid](jcb+j,i);
        bx = U->Q[bxid](jcb+j,i);
        by = U->Q[byid](jcb+j,i);

        U->Q[uid](jgb,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jgb,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jgb,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jgb,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jgb,i) = 2.0*U->Q[rhoid](jgb+1,i) - U->Q[rhoid](jgb+2,i);
        U->Q[wid](jgb,i) = 2.0*U->Q[wid](jgb+1,i) - U->Q[wid](jgb+2,i);
        U->Q[pid](jgb,i) = 2.0*U->Q[pid](jgb+1,i) - U->Q[pid](jgb+2,i);
        U->Q[bzid](jgb,i) = 2.0*U->Q[bzid](jgb+1,i) - U->Q[bzid](jgb+2,i);
      }
    // Is no need to send data because bottom is physical boundary
    //printf("My bottom on phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }
  else // not physical boundary send my stuff on bottom to the rank below me
  {
    int jcb = C.num_ghost; // first cell layer on bottom
    for(int j = 0; j < C.num_ghost; j++)
    {
      // no flip sig
      tempBT->Q[rhoid].row(j) = U->Q[rhoid].row(jcb+j);
      tempBT->Q[wid].row(j) = U->Q[wid].row(jcb+j);
      tempBT->Q[bzid].row(j) = U->Q[bzid].row(jcb+j);
      tempBT->Q[pid].row(j) = U->Q[pid].row(jcb+j);
      tempBT->Q[uid].row(j) =  U->Q[uid].row(jcb+j);
      tempBT->Q[vid].row(j) =  U->Q[vid].row(jcb+j);
      tempBT->Q[bxid].row(j) = U->Q[bxid].row(jcb+j);
      tempBT->Q[byid].row(j) = U->Q[byid].row(jcb+j);
    }
    // take my data and send it to the processor below me which is not phys boundary
    MPI_Send(tempBT->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, d, 444, com2d); // tag 444 is top
    //printf("My bottom no phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }

  if(u < 0) // physical boundary on top no need send data
  {
    int jft = njx.rows()-1; // top face
    int jgt;
    int iff;

    int jct;
    double uvel, vvel, bx, by, nx,ny;

    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {
        iff = i-C.num_ghost;
        jct = njc - C.num_ghost - 1 - j;
        jgt = (njc-C.num_ghost)+j;

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

      }
    // I no need send data directly modify u with phys bc
    //printf("My top on phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }
  else
  {
    int jct = (njc-C.num_ghost)-1; // first cell layer on top
    for(int j = 0; j < C.num_ghost; j++)
    {
      // no flip sig
      tempBT->Q[rhoid].row(j) = U->Q[rhoid].row(jct-j);
      tempBT->Q[wid].row(j) = U->Q[wid].row(jct-j);
      tempBT->Q[bzid].row(j) = U->Q[bzid].row(jct-j);
      tempBT->Q[pid].row(j) = U->Q[pid].row(jct-j);
      tempBT->Q[uid].row(j) =  U->Q[uid].row(jct-j);
      tempBT->Q[vid].row(j) =  U->Q[vid].row(jct-j);
      tempBT->Q[bxid].row(j) = U->Q[bxid].row(jct-j);
      tempBT->Q[byid].row(j) = U->Q[byid].row(jct-j);
    }
    // take my data and send it to the processor above me which is not phys boundary
    MPI_Send(tempBT->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, u, 333, com2d); // tag 333 is bottom
    //printf("My bottom no phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);

    //printf("My top no phys my rank is %d: My neighbors are left:%d right:%d up:%d down:%d\n", rank, l, r, u ,d);
  }
  delete[] tempBT->Q_raw; tempBT->Q_raw = NULL;
}

void setBcRecv(Map2Eigen* U, MPI_Comm& com2d,  constants C)
{
  int rank;
  int size;
  MPI_Comm_rank(com2d, &rank);
  MPI_Comm_size(com2d, &size);

  int u, l, d, r;
  computeCoord(l, r, d, u, com2d);

  int njc = U->Q[rhoid].rows();
  int nic = U->Q[rhoid].cols();

  MPI_Request request; 
  MPI_Status status; 

  // lets receive left and right if we need to
  ///////////////////// LEFT RIGHT ////////////////////////
  Map2Eigen *tempL= new Map2Eigen(njc , C.num_ghost, NEQ);
  Map2Eigen *tempR= new Map2Eigen(njc , C.num_ghost, NEQ);

  // Receive data for the left
  MPI_Recv(tempL->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, l, 111, com2d, &status);
  MPI_Recv(tempR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, r, 222, com2d, &status);

  // Fill in data for the left
  for(int eq = 0; eq < NEQ; eq++)
  {
    int ig=C.num_ghost-1;
    for(int i = 0; i < C.num_ghost; i++)
      U->Q[eq].col(ig-i) = tempL->Q[eq].col(i);

    ig = (nic-C.num_ghost);
    for(int i = 0; i < C.num_ghost; i++)
      U->Q[eq].col(ig+i) = tempR->Q[eq].col(i);
  }
  delete[] tempL->Q_raw, tempL->Q_raw=NULL;
  delete[] tempR->Q_raw, tempR->Q_raw=NULL;


  ///////////////////// TOP BOTTOM ////////////////////////
  Map2Eigen *tempB= new Map2Eigen(C.num_ghost, nic, NEQ);
  Map2Eigen *tempT= new Map2Eigen(C.num_ghost, nic, NEQ);


  if (d > 0)// then need receive lower ghost info
  {
    MPI_Recv(tempB->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, d, 333, com2d, &status);
    int jg = C.num_ghost-1;
    for(int eq = 0; eq < NEQ; eq++)
      for(int j = 0; j < C.num_ghost; j++)
        U->Q[eq].row(jg-j) = tempB->Q[eq].row(j);
  }

  if (u > 0)// then need receive
  {

    MPI_Recv(tempT->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, u, 444, com2d, &status);

    int jg = (njc-C.num_ghost);
    for(int eq = 0; eq < NEQ; eq++)
      for(int j = 0; j < C.num_ghost; j++)
        U->Q[eq].row(jg+j) = tempT->Q[eq].row(j);
  }
  delete[] tempB->Q_raw; tempB->Q_raw=NULL;
  delete[] tempT->Q_raw; tempT->Q_raw=NULL;

  MPI_Barrier(com2d);
}

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

  int ic, sign, out, tag; // create items to be passed to if statements
  int sendSize = njc*C.num_ghost*NEQ;
  /////////////////// do left////////////////
  ic = C.num_ghost;
  sign = +1;
  out = l;
  tag = 222; // right recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempLR->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(l>rank) // phys bndry else not
  {
    // flip sign
    tempLR->Q[uid]  =-1.0*tempLR->Q[uid];
    tempLR->Q[vid]  =-1.0*tempLR->Q[vid];
    tempLR->Q[bxid] =-1.0*tempLR->Q[bxid];
    tempLR->Q[byid] =-1.0*tempLR->Q[byid];
  }
  MPI_Isend(tempLR->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request);//tag222 is right
  ////////////////////////////////////////

  ////////////// do right ////////////////
  ic =(nic-C.num_ghost)-1; // first cell layer on right 
  sign = -1;
  out = r;
  tag = 111; // left recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempLR->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(r < rank)// physical boundary on right  need flip
  {
    tempLR->Q[uid]  =-1.0*tempLR->Q[uid];
    tempLR->Q[vid]  =-1.0*tempLR->Q[vid];
    tempLR->Q[bxid] =-1.0*tempLR->Q[bxid];
    tempLR->Q[byid] =-1.0*tempLR->Q[byid];
  }
  MPI_Isend(tempLR->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag222 is right
  /////////////////////////////////////
  delete[] tempLR->Q_raw; tempLR->Q_raw = NULL;

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// SEND BOTT AND TOP  DATA  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  Map2Eigen *tempBT  = new Map2Eigen(C.num_ghost, nic, NEQ);
  sendSize = nic*C.num_ghost*NEQ;

  /////////////////// do bottom ///////////////////////
  int jfo, jgo, jco;
  jfo = 0;// bottom face
  jco = C.num_ghost; // bott inter offset
  jgo = C.num_ghost- 1;
  sign = +1;

  int iff;
  int jc;
  int jg;
  double uvel, vvel, bx, by, nx,ny;

  if(d < 0) // phys bndry
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {
        iff = i - C.num_ghost;// index for normvecs
        jc = jco + sign*j; // cell interior
        jg = jgo - sign*j; // ghost cell row-layer

        // bottom boundary no penetration conducting wall
        nx = njx(jfo,iff);
        ny = njy(jfo,iff);
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);

        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[pid](jg,i) =   2.0*U->Q[pid](jg+sign*1,i) - U->Q[pid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry
  {
    out = d;
    tag = 444;
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempBT->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    MPI_Isend(tempBT->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag444top
  }


  ////////////////// do upper //////////////////////
  jfo = njx.rows()-1;// top face
  jco = njc - C.num_ghost - 1; // top inter offset
  jgo = (njc-C.num_ghost);
  sign = -1;
  out = u;
  tag = 333;

  if (u < 0) // phys bndry on top
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {
        iff = i - C.num_ghost;// index for normvecs
        jc = jco + sign*j; // cell interior
        jg = jgo - sign*j; // ghost cell row-layer

        // bottom boundary no penetration conducting wall
        nx = njx(jfo,iff);
        ny = njy(jfo,iff);
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);

        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[pid](jg,i) =   2.0*U->Q[pid](jg+sign*1,i) - U->Q[pid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry on top
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempBT->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    MPI_Isend(tempBT->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag333bot
  }
  //////////////////////////////////////////////

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
  int recvSize = njc*C.num_ghost*NEQ;

  // Receive data for the left
  //MPI_Recv(tempL->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, l, 111, com2d, &status);
  //MPI_Recv(tempR->Q_raw, njc*C.num_ghost*NEQ, MPI_DOUBLE, r, 222, com2d, &status);
  MPI_Irecv(tempL->Q_raw, recvSize, MPI_DOUBLE, l, 111, com2d, &request);
  //MPI_Wait(&request, &status);
  MPI_Irecv(tempR->Q_raw, recvSize, MPI_DOUBLE, r, 222, com2d, &request);
  //MPI_Wait(&request, &status);
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

  recvSize = nic*C.num_ghost*NEQ;

  if (d > 0)// then need receive lower ghost info
  {
    //MPI_Recv(tempB->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, d, 333, com2d, &status);
    MPI_Irecv(tempB->Q_raw, recvSize, MPI_DOUBLE, d, 333, com2d, &request);
    MPI_Wait(&request, &status);
    int jg = C.num_ghost-1;
    for(int eq = 0; eq < NEQ; eq++)
      for(int j = 0; j < C.num_ghost; j++)
        U->Q[eq].row(jg-j) = tempB->Q[eq].row(j);
  }

  if (u > 0)// then need receive
  {

    //MPI_Recv(tempT->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, u, 444, com2d, &status);
    MPI_Irecv(tempT->Q_raw, nic*C.num_ghost*NEQ, MPI_DOUBLE, u, 444, com2d, &request);

    MPI_Wait(&request, &status);
    int jg = (njc-C.num_ghost);
    for(int eq = 0; eq < NEQ; eq++)
      for(int j = 0; j < C.num_ghost; j++)
        U->Q[eq].row(jg+j) = tempT->Q[eq].row(j);
  }
  delete[] tempB->Q_raw; tempB->Q_raw=NULL;
  delete[] tempT->Q_raw; tempT->Q_raw=NULL;

  //MPI_Barrier(com2d);
}




void mpiSetBc(Map2Eigen* U, const RowMajorMatrixXd& nix, const RowMajorMatrixXd& niy, const RowMajorMatrixXd& njx, const RowMajorMatrixXd& njy, MPI_Comm& com2d,  constants C)
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
  MPI_Request requestIn[4]; 
  MPI_Status status; 

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// SEND LEFT AND RIGHT DATA  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  Map2Eigen *tempLR  = new Map2Eigen(njc , C.num_ghost, NEQ);

  int ic, sign, out, tag; // create items to be passed to if statements
  int sendSize = njc*C.num_ghost*NEQ;
  /////////////////// do left////////////////
  ic = C.num_ghost;
  sign = +1;
  out = l;
  tag = 222; // right recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempLR->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(l>rank) // phys bndry else not
  {
    // flip sign
    tempLR->Q[uid]  =-1.0*tempLR->Q[uid];
    tempLR->Q[vid]  =-1.0*tempLR->Q[vid];
    tempLR->Q[bxid] =-1.0*tempLR->Q[bxid];
    tempLR->Q[byid] =-1.0*tempLR->Q[byid];
  }
  // sent to left 
  MPI_Isend(tempLR->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request);//tag222 is right

  // recv from left fill U
  int in = l;
  tag = 111;
  MPI_Irecv(tempLR->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn[0]);

  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      U->Q[eq].col(ic+sign*i) = tempLR->Q[eq].col(i);
  ////////////////////////////////////////

  ////////////// do right ////////////////
  ic =(nic-C.num_ghost)-1; // first cell layer on right 
  sign = -1;
  out = r;
  tag = 111; // left recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempLR->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(r < rank)// physical boundary on right  need flip
  {
    tempLR->Q[uid]  =-1.0*tempLR->Q[uid];
    tempLR->Q[vid]  =-1.0*tempLR->Q[vid];
    tempLR->Q[bxid] =-1.0*tempLR->Q[bxid];
    tempLR->Q[byid] =-1.0*tempLR->Q[byid];
  }
  // sent to right
  MPI_Isend(tempLR->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag222 is right
  
  tag = 222;
  in = r;
  // recv from right
  MPI_Irecv(tempLR->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn[1]);
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      U->Q[eq].col(ic+sign*i) = tempLR->Q[eq].col(i);

  /////////////////////////////////////
  delete[] tempLR->Q_raw; tempLR->Q_raw = NULL;

  /////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// SEND BOTT AND TOP  DATA  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////
  Map2Eigen *tempBT  = new Map2Eigen(C.num_ghost, nic, NEQ);
  sendSize = nic*C.num_ghost*NEQ;

  /////////////////// do bottom ///////////////////////
  int jfo, jgo, jco;
  jfo = 0;// bottom face
  jco = C.num_ghost; // bott inter offset
  jgo = C.num_ghost- 1;
  sign = +1;

  int iff;
  int jc;
  int jg;
  double uvel, vvel, bx, by, nx,ny;

  if(d < 0) // phys bndry
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {
        iff = i - C.num_ghost;// index for normvecs
        jc = jco + sign*j; // cell interior
        jg = jgo - sign*j; // ghost cell row-layer

        // bottom boundary no penetration conducting wall
        nx = njx(jfo,iff);
        ny = njy(jfo,iff);
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);

        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[pid](jg,i) =   2.0*U->Q[pid](jg+sign*1,i) - U->Q[pid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry
  {
    out = d;
    tag = 444;
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempBT->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    // send bott dat
    MPI_Isend(tempBT->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag444top

    in = d;
    tag = 333;
    // recv bott dat
    MPI_Irecv(tempBT->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn[3]); //tag444top

    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q[eq].row(jco+sign*j) = tempBT->Q[eq].row(j);
  }

  ////////////////// do upper //////////////////////
  jfo = njx.rows()-1;// top face
  jco = njc - C.num_ghost - 1; // top inter offset
  jgo = (njc-C.num_ghost);
  sign = -1;
  out = u;
  tag = 333;

  if (u < 0) // phys bndry on top
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int i = C.num_ghost; i < nic-C.num_ghost; i++)
      {
        iff = i - C.num_ghost;// index for normvecs
        jc = jco + sign*j; // cell interior
        jg = jgo - sign*j; // ghost cell row-layer

        // bottom boundary no penetration conducting wall
        nx = njx(jfo,iff);
        ny = njy(jfo,iff);
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);

        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[pid](jg,i) =   2.0*U->Q[pid](jg+sign*1,i) - U->Q[pid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry on top
  {
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempBT->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    // send top data
    MPI_Isend(tempBT->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &request); //tag333bot

    in = u;
    tag = 444;
    // recv top data
    
    MPI_Irecv(tempBT->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn[4]); //tag444top

    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q[eq].row(jco+sign*j) = tempBT->Q[eq].row(j);
  }
  //////////////////////////////////////////////

  delete[] tempBT->Q_raw; tempBT->Q_raw = NULL;
}


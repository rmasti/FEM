#include "mhdRT.hpp"
 
void mpiSetBc(
    Map2Eigen* U,                 //Out: conserved vars
    const RowMajorMatrixXd& nix,  //In: norm vec xcomp i dir
    const RowMajorMatrixXd& niy,  //In: norm vec ycomp i dir
    const RowMajorMatrixXd& njx,  //In: norm vec xcomp j dir
    const RowMajorMatrixXd& njy,  //In: norm vec ycomp j dir
    MPI_Comm& com2d,              //In: Communicator 
    constants C)                   
{
  int rank;
  int size;
  MPI_Comm_rank(com2d, &rank);
  MPI_Comm_size(com2d, &size);

  int u, l, d, r;
  computeCoord(l, r, d, u, com2d);

  int njc = U->Q[rhoid].rows();
  int nic = U->Q[rhoid].cols();

  MPI_Request requestOut[4]; 
  MPI_Request requestIn; 
  MPI_Status status; 

  //////////////////////////// SEND LEFT AND RIGHT DATA  //////////////////////////////////

  int ic, sign, out, tag, in; // create items to be passed to if statements
  int sendSize = njc*C.num_ghost*NEQ;
  /////////////////// do left////////////////
  Map2Eigen *tempLout = new Map2Eigen(njc , C.num_ghost, NEQ);
  ic = C.num_ghost;
  sign = +1;
  out = l;
  tag = 222; // right recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempLout->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(l>rank) // phys bndry else not
  {
    // flip sign
    tempLout->Q[uid]  =-1.0*tempLout->Q[uid];
    tempLout->Q[vid]  =-1.0*tempLout->Q[vid];
    tempLout->Q[bxid] =-1.0*tempLout->Q[bxid];
    tempLout->Q[byid] =-1.0*tempLout->Q[byid];
    //printf("My rank is %d; l = %d, r = %d, d = %d, u = %d\n", rank, l, r, d, u);
  }
  // sent to left 
  MPI_Isend(tempLout->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &requestOut[0]);//tag222 is right

  ////////////// do right ////////////////
  Map2Eigen *tempRout = new Map2Eigen(njc , C.num_ghost, NEQ);
  ic =(nic-C.num_ghost)-1; // first cell layer on right 
  sign = -1;
  out = r;
  tag = 111; // left recv
  // fill lr
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      tempRout->Q[eq].col(i) = U->Q[eq].col(ic+sign*i);

  if(r < rank)// physical boundary on right  need flip
  {
    tempRout->Q[uid]  =-1.0*tempRout->Q[uid];
    tempRout->Q[vid]  =-1.0*tempRout->Q[vid];
    tempRout->Q[bxid] =-1.0*tempRout->Q[bxid];
    tempRout->Q[byid] =-1.0*tempRout->Q[byid];
  }
  // sent to right
  MPI_Isend(tempRout->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &requestOut[1]); //tag222 is right

  //////////////////////////// SEND BOTT AND TOP  DATA  //////////////////////////////////
  sendSize = nic*C.num_ghost*NEQ;

  /////////////////// do bottom ///////////////////////
  Map2Eigen *tempBout  = new Map2Eigen(C.num_ghost, nic, NEQ);
  Map2Eigen *tempBin  = new Map2Eigen(C.num_ghost, nic, NEQ);
  int jfo, jgo, jco;
  jfo = 0;// bottom face
  jco = C.num_ghost; // bott inter offset
  jgo = C.num_ghost - 1;
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
        /*
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);
        */
        uvel = U->Q[uid](jc,i);
        vvel = U->Q[vid](jc,i);
        bx =  U->Q[bxid](jc,i);
        by =  U->Q[byid](jc,i);
 
        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate from cells above
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[piid](jg,i) =   2.0*U->Q[piid](jg+sign*1,i) - U->Q[piid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry
  {

    // send bott dat
    out = d;
    tag = 444;
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempBout->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    MPI_Isend(tempBout->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &requestOut[2]); //tag444top
  }

  ////////////////// do upper //////////////////////
  Map2Eigen *tempTout  = new Map2Eigen(C.num_ghost, nic, NEQ);
  Map2Eigen *tempTin  = new Map2Eigen(C.num_ghost, nic, NEQ);
  jfo = njx.rows()-1;// top face
  jco = njc - C.num_ghost - 1; // top inter offset
  jgo = (njc-C.num_ghost);
  sign = -1;

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

        /**** BELIEVED BUG *****
        uvel = U->Q[uid](jc+sign*j,i);
        vvel = U->Q[vid](jc+sign*j,i);
        bx =  U->Q[bxid](jc+sign*j,i);
        by =  U->Q[byid](jc+sign*j,i);
        */
        uvel = U->Q[uid](jc,i);
        vvel = U->Q[vid](jc,i);
        bx =  U->Q[bxid](jc,i);
        by =  U->Q[byid](jc,i);

        U->Q[uid](jg,i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        U->Q[vid](jg,i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        U->Q[bxid](jg,i) = -nx*(bx*nx+by*ny) - ny*(-bx*ny + by*nx);
        U->Q[byid](jg,i) = -ny*(bx*nx+by*ny) + nx*(-bx*ny + by*nx);
        // extrapolate
        U->Q[rhoid](jg,i) = 2.0*U->Q[rhoid](jg+sign*1,i) - U->Q[rhoid](jg+sign*2,i);
        U->Q[wid](jg,i) =   2.0*U->Q[wid](jg+sign*1,i) - U->Q[wid](jg+sign*2,i);
        U->Q[piid](jg,i) =   2.0*U->Q[piid](jg+sign*1,i) - U->Q[piid](jg+sign*2,i);
        U->Q[bzid](jg,i) =  2.0*U->Q[bzid](jg+sign*1,i) - U->Q[bzid](jg+sign*2,i);
      }
  }
  else // not phys bndry on top
  {
    // send top data
    out = u;
    tag = 333;
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        tempTout->Q[eq].row(j) = U->Q[eq].row(jco+sign*j);

    MPI_Isend(tempTout->Q_raw, sendSize, MPI_DOUBLE, out, tag, com2d, &requestOut[3]); //tag333bot

  }

  /////////////////////// RECEIVE DATA ///////////////////////  
  // recv from left fill U
  Map2Eigen *tempLin  = new Map2Eigen(njc , C.num_ghost, NEQ);
  sendSize = njc*C.num_ghost*NEQ;
  in = l;
  tag = 111;
  ic = C.num_ghost;
  sign = +1;
  MPI_Irecv(tempLin->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn);
  MPI_Wait(&requestIn, &status);
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      U->Q[eq].col((ic-1)-sign*i) = tempLin->Q[eq].col(i);
  ////////////////////////////////////////
  // recv from right
  Map2Eigen *tempRin  = new Map2Eigen(njc , C.num_ghost, NEQ);
  sendSize = njc*C.num_ghost*NEQ;
  in = r;
  tag = 222;
  ic =(nic-C.num_ghost)-1; // first cell layer on right 
  sign = -1;
  MPI_Irecv(tempRin->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn);
  MPI_Wait(&requestIn, &status);
  for(int i = 0; i < C.num_ghost; i++)
    for(int eq = 0; eq < NEQ; eq++)
      U->Q[eq].col(ic+1-sign*i) = tempRin->Q[eq].col(i); // fill in data on right
  /////////////////////////////////////
  // recv bott dat
  if (d >= 0)
  {
    Map2Eigen *tempBin  = new Map2Eigen(C.num_ghost, nic, NEQ);
    sendSize = C.num_ghost*nic*NEQ;
    in = d;
    tag = 333;
    jco = C.num_ghost; // bott inter offset
    sign = +1;
    MPI_Irecv(tempBin->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn); //tag444top
    MPI_Wait(&requestIn, &status);
    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q[eq].row(jco-1-sign*j) = tempBin->Q[eq].row(j); // fill in bottom data shift and flip sign
    delete tempBin; tempBin = NULL;
  }
  /////////////////////////////////////
  // recv top data
  if(u >= 0)
  {
    Map2Eigen *tempTin  = new Map2Eigen(C.num_ghost, nic, NEQ);
    sendSize = C.num_ghost*nic*NEQ;
    in = u;
    tag = 444; // data coming in from top
    jco = njc - C.num_ghost - 1; // top inter offset
    sign = -1;
    MPI_Irecv(tempTin->Q_raw, sendSize, MPI_DOUBLE, in, tag, com2d, &requestIn); //tag444top
    MPI_Wait(&requestIn, &status);

    for(int j = 0; j < C.num_ghost; j++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q[eq].row(jco+1-sign*j) = tempTin->Q[eq].row(j);
    delete tempTin; tempTin= NULL;
  }
  ///////////////////////////////////

  delete tempTin; tempTin = NULL;
  delete tempTout; tempTout = NULL;

  //MPI_Barrier(com2d);
  delete tempLin; tempLin = NULL;
  delete tempLout; tempLout = NULL;

  delete tempRin; tempRin = NULL;
  delete tempRout; tempRout = NULL;

  delete tempBin; tempBin = NULL;
  delete tempBout; tempBout = NULL;
}

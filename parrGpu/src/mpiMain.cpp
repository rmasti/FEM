#include "mhdRT.hpp"
#include "occa.hpp"
//#include <cuda_runtime_api.h>

int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv); 

  constants C;
  C.f_limiter = 5;
  C.num_ghost = 3;
  C.cfl = 0.45;
  C.nmax = 1000000;
  C.wint = 1000;
  C.pint = 100;
  double A = (2-1)/(1.0+2.0);
  double tend = 6.0/sqrt(A*ACCEL*2);
  double dt;
  //occa::umalloc(sizeof(double), &dt);
  MatrixXd time(1, C.nmax);
  time(0,0)=0.0;

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  string mesh = "../mesh/560x480.msh"; //debugMatlab.msh";
  //string mesh = "../mesh/360x300.msh"; //debugMatlab.msh";
  //string mesh = "../mesh/16x10.msh"; //debugMatlab.msh";
  string outputFolder = "./output/";
  //string outputFolder = "/mnt/c/Users/rlm78/Downloads/FEM/";

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

  int coordMax[2]={1, 1};
  MPI_Cart_coords(com2d, rank, 2, coordMax);

  MPI_Allreduce(MPI_IN_PLACE, &coordMax[0], 1, MPI_INT, MPI_MAX, com2d);
  MPI_Allreduce(MPI_IN_PLACE, &coordMax[1], 1, MPI_INT, MPI_MAX, com2d);

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
  primToCons(U, V);

  U_RK = U;
  outputArrayMap(outputFolder, "UL", U, rank);
 
  // send and recv func
  //setBcSend(U, nixL, niyL, njxL, njyL, com2d,  C);
  //setBcRecv(U, com2d, C);
  //MPI_Barrier(com2d);

  // recv while you send
  //cout << " setting bc " << endl;
   
  
  if(rank == 0)
    cout << " Setting BC " << endl;


  mpiSetBc(U, nixL, niyL, njxL, njyL, com2d,  C);

  if(rank == 0)
    cout << " Set BC " << endl;
 
  MPI_Barrier(com2d);

  //outputArrayMap(outputFolder, "UL", U, rank);
  //outputArrayMap(outputFolder, "VL", V, rank);
  outputArray(outputFolder, "xcL_g", xcL_g, rank+size);

  dt = computeTimeStepU(VolumeL, AiL, AjL, nixL, niyL, njxL, njyL, U, C);
  MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, com2d);

  MPI_Barrier(com2d);

  ////////// load on to device ////////// 
  // Things needed to put on to device
  // nixL, niyL 
  // njxL, njyL
  // AjL, AiL
  // Volume L
  // U, U_RK
  // U_B, U_T, U_L, U_R, F, G
  // S, Res, xcL, ycL
  
  // total bytes
  
  int bytes = njc_g*nic_g*NEQ*2; // U, U_RK
  bytes += njc*(nic+1)*(NEQ*3+3); // U_B, U_T, G, nixL, niyL, AiL
  bytes += (njc+1)*(nic)*(NEQ*3+3); // U_L, U_R, F, njxL, njyL, AjL
  bytes += njc*nic*(3+NEQ*2); // vol, xcL, ycL, S, Res
  bytes = bytes*8; 
  
  if(rank==0)
    printf("total gb loaded onto gpu: tot = %lf\n", bytes*1.0e-9);

  
  occa::device device;
  char setupString[BUFSIZ];
  sprintf(setupString, "mode: 'CUDA', device_id : %d", 0); // odd ranks to 1 even ranks to 0
  //sprintf(setupString, "mode: 'Serial', device_id : %d", 0); // odd ranks to 1 even ranks to 0
  device.setup(setupString); 

  occa::memory o_U = device.malloc(njc_g*nic_g*NEQ*sizeof(double), U->Q_raw);
  occa::memory o_U_RK = device.malloc(njc_g*nic_g*NEQ*sizeof(double), U_RK->Q_raw);

  occa::memory o_njxL = device.malloc((nic)*(njc+1)*sizeof(double), njxL.data());
  occa::memory o_njyL = device.malloc((nic)*(njc+1)*sizeof(double), njyL.data());
  occa::memory o_AjL =  device.malloc((nic)*(njc+1)*sizeof(double), AjL.data());
  occa::memory o_U_B =  device.malloc(nic*(njc+1)*NEQ*sizeof(double), U_B->Q_raw);
  occa::memory o_U_T =  device.malloc(nic*(njc+1)*NEQ*sizeof(double), U_T->Q_raw);
  occa::memory o_G =    device.malloc(nic*(njc+1)*NEQ*sizeof(double), G->Q_raw);

  occa::memory o_nixL = device.malloc(njc*(nic+1)*sizeof(double), nixL.data());
  occa::memory o_niyL = device.malloc(njc*(nic+1)*sizeof(double), niyL.data());
  occa::memory o_AiL =  device.malloc(njc*(nic+1)*sizeof(double), AiL.data()); 
  occa::memory o_U_L =  device.malloc((nic+1)*(njc)*NEQ*sizeof(double), U_L->Q_raw);
  occa::memory o_U_R =  device.malloc((nic+1)*(njc)*NEQ*sizeof(double), U_R->Q_raw);
  occa::memory o_F =    device.malloc((nic+1)*(njc)*NEQ*sizeof(double), F->Q_raw);
  
  occa::memory o_xcL = device.malloc((njc)*(nic)*sizeof(double), xcL.data());
  occa::memory o_ycL = device.malloc((njc)*(nic)*sizeof(double), ycL.data());
  occa::memory o_VolumeL = device.malloc((njc)*(nic)*sizeof(double), VolumeL.data());
  occa::memory o_S = device.malloc((njc)*(nic)*NEQ*sizeof(double), S->Q_raw);
  occa::memory o_Res = device.malloc((njc)*(nic)*NEQ*sizeof(double), Res->Q_raw);


  // memory needed for rk staging

  occa::properties props;
  props["defines/o_NEQ"]=NEQ;
  props["defines/o_ACCEL"]=ACCEL;
  props["defines/o_ng"]=C.num_ghost;
  props["defines/o_rhoid"]=rhoid;
  props["defines/o_uid"]=uid;
  props["defines/o_vid"]=vid;
  props["defines/o_wid"]=wid;
  props["defines/o_piid"]=piid;
  props["defines/o_bxid"]=bxid;
  props["defines/o_byid"]=byid;
  props["defines/o_bzid"]=bzid;
  props["defines/o_njc_g"]=njc_g;
  props["defines/o_nic_g"]=nic_g;
  props["defines/o_limiter"]=C.f_limiter;
  props["defines/o_gamma"]=GAMMA;
  props["defines/o_mu0"]=MU0;
  props["defines/o_rkorder"]=RKORDER;
  props["defines/o_njc"]=njc;
  props["defines/o_nic"]=nic;

  //props["includes/"]+="include/muscl.hpp";

  //occa::kernel occaMap2eigen = device.buildKernel("test/occaMap2Eigen.okl", "occaMap2Eigen");//, props);
  
  cout << " Building Kernels" << endl;
  occa::kernel computeSourceTerm = device.buildKernel("src/computeSourceTerm.okl", "computeSourceTerm", props);
  occa::kernel MUSCL_LR = device.buildKernel("src/MUSCL.okl", "MUSCL_LR", props);
  occa::kernel MUSCL_BT = device.buildKernel("src/MUSCL.okl", "MUSCL_BT", props);
  occa::kernel compute2dFluxF = device.buildKernel("src/compute2dFlux.okl", "compute2dFluxF", props);
  occa::kernel compute2dFluxG = device.buildKernel("src/compute2dFlux.okl", "compute2dFluxG", props);
  occa::kernel computeRes = device.buildKernel("src/computeRes.okl", "computeRes", props);
  occa::kernel rungeKutta = device.buildKernel("src/rungeKutta.okl", "rungeKutta", props);

  device.finish();

  int n=0;

  outputArrayMap(outputFolder, "UL", U, rank);
  stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);

  clock_t t;
  float avgT[6];
  if(rank == 0)
    cout << " Entering Time Loop " << endl;

  U_RK = U;
  while(time(0, n) < tend)// && n < 10)
  {
    for (int k = 0; k < RKORDER; k++)
    {
      t = clock();
      //cout << "here" << endl;
      mpiSetBc(U, nixL, niyL, njxL, njyL, com2d,  C);
      t = clock()-t;
      avgT[0] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);
      MPI_Barrier(com2d);
      o_U_RK.copyFrom(U_RK->Q_raw);
      device.finish();

      t = clock();
      //computeSourceTerm(S, U_RK, xcL, ycL, C);
      computeSourceTerm(o_S, o_U_RK, o_xcL, o_ycL);
      device.finish();
      //o_S.copyTo(S->Q_raw);
      //device.finish();
      //device.finish();
      //cout << S->Q[uid] << endl; t = clock()-t;
      t = clock()-t;
      avgT[1] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);

      t = clock();
      MUSCL_LR(o_U_L, o_U_R, o_U_RK);
      device.finish();
      MUSCL_BT(o_U_B, o_U_T, o_U_RK);
      device.finish();
      //o_U_L.copyTo(U_L->Q_raw);
      //device.finish();
      //o_U_R.copyTo(U_R->Q_raw);
      //device.finish();
      //o_U_B.copyTo(U_B->Q_raw);
      //device.finish();
      //o_U_T.copyTo(U_T->Q_raw);
      //device.finish();
      //cout << endl;
      //cout << U->Q[rhoid] << endl; 
      //cout << U_R->Q[rhoid] << endl; 
      //cout << U_T->Q[rhoid] << endl; 
      //cout << U_B->Q[rhoid] << endl; 
      //cout << U_T->Q[rhoid] << endl; 

      t = clock()-t;
      avgT[2] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);

      t = clock();
      compute2dFluxF(o_F, o_U_L, o_U_R, o_nixL, o_niyL);
      device.finish();
      compute2dFluxG(o_G, o_U_B, o_U_T, o_njxL, o_njyL);
      device.finish();

      //o_F.copyTo(F->Q_raw);
      //device.finish();
      //o_G.copyTo(G->Q_raw);

      //cout << G->Q[rhoid] << endl;
      //cout << G->Q[rhoid] << endl;

      t = clock()-t;
      avgT[3] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);

      t = clock();
      //computeRes(Res, S, F, G, AjL, AiL, VolumeL, C);
      computeRes(o_Res, o_S, o_F, o_G, o_AjL, o_AiL, o_VolumeL);
      device.finish();
      //o_Res.copyTo(Res->Q_raw);
      //device.finish();
      //cout << Res->Q[rhoid] << endl;
      

      t = clock()-t;
      avgT[4] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);

      t = clock();
      //rungeKutta(U_RK, U, Res, VolumeL, k, dt, C);
      //o_dt.copyFrom(*dt);
      //o_k.copyFrom(*k);
      rungeKutta(o_U_RK, o_U, o_Res, o_VolumeL, k, dt);
      device.finish();
      o_U_RK.copyTo(U_RK->Q_raw);
      device.finish();

      //cout << U_RK->Q[rhoid] << endl;

      t = clock()-t;
      avgT[5] =  ((float)t)/(1.0e-3*CLOCKS_PER_SEC);
    }
    //cout << " Made it past rk stages " << endl;
  
    U=U_RK;
    //cout << U->Q[rhoid]<< endl;

    if(time(0,n)+dt > tend)// end at exact time
      dt = tend-time(0,n);
    else
    {
      dt = computeTimeStepU(VolumeL, AiL, AjL, nixL, niyL, njxL, njyL, U, C);
      MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MIN, com2d);
    }

    time(0,n+1) = time(0,n) + dt;
    n++;

    if(n%C.pint == 0)
      if (rank == 0){
        cout << "time  = " << time(0,n) << ", n = " << n << endl;
        printf("setBc=%lf, srcTerm=%lf, muscl=%lf, 2dflux=%lf, res=%lf, rk=%lf\n", avgT[0],avgT[1],avgT[2],avgT[3],avgT[4],avgT[5]);
      }

    if(n%C.wint == 0) // Write to file by rank 0 gather and output
    {
      if(rank == 0){
        cout << " Write To File..." << endl;
      }

      stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);
    }
  }
  stitchMap2EigenWrite(outputFolder, "U", U, n,coordMax, com2d, C);
  //occaFree(MUSCL_LR);
  //occaFree(MUSCL_BT);
  //occaFree(o_U);
  //occaFree(props);
  //outputArrayMap(outputFolder, "UL", U, rank);

  delete V; V = NULL;
  delete U; U = NULL;
  //delete U_RK; U_RK= NULL;
  delete Res; Res= NULL;
  delete S; S= NULL;
  delete F; F= NULL;
  delete G; G= NULL;
  delete U_B; U_B= NULL;
  delete U_T; U_T= NULL;
  delete U_L; U_L= NULL;
  delete U_R; U_R= NULL;

  occa::free(device);

  MPI_Finalize();
  return 0;

}

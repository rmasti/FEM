#include "mhdRT.hpp"

MPI_Comm meshBlock(const string mesh, const string outputFolder, RowMajorMatrixXd &xcL_g, RowMajorMatrixXd &ycL_g, RowMajorMatrixXd &nixL, RowMajorMatrixXd &niyL, RowMajorMatrixXd &njxL, RowMajorMatrixXd &njyL, RowMajorMatrixXd &AiL, RowMajorMatrixXd &AjL, RowMajorMatrixXd &VolumeL, constants C)
{
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm com2d;

  int ni, nj;

  int cnt=size;
  int n=0;

  while(cnt >=2)
  {
    cnt = cnt/2;
    n++;
  }
  nj = pow(2, (n/2)); // processors
  ni = pow(2, (n/2) + n%2);

  int dim[2] = {nj, ni};
  int reorder = TRUE;//FALSE;
  int period[2];
  period[1] = TRUE;
  period[0] = FALSE;
  int coord[2];
  int id;

  int u, d, l, r;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &com2d);
  MPI_Comm_rank(com2d, &rank);
  //MPI_Cart_coords(com2d, rank, 2, coord);
  //  printf("My rank is %d: My coordinates are %d, %d\n", rank, coord[0], coord[1]);

  int nb[2];
  // All the variables geometry wise needed locally to loop
  /*
  RowMajorMatrixXd xcL_g, ycL_g;
  RowMajorMatrixXd nixL, niyL;
  RowMajorMatrixXd njxL, njyL;
  RowMajorMatrixXd AjL, AiL;
  RowMajorMatrixXd VolumeL;
  */


  // have the first processor read in data and give it to all other ranks
  if(rank == 0)
  {
    RowMajorMatrixXd xn;
    RowMajorMatrixXd yn;
    RowMajorMatrixXd xc, yc;
    // Read in mesh 
    inputMesh(xn,yn,xc,yc, mesh);

    // extrap to ghost cells
    RowMajorMatrixXd xc_g(xc.rows()+2*C.num_ghost, xc.cols()+2*C.num_ghost);
    RowMajorMatrixXd yc_g(yc.rows()+2*C.num_ghost, yc.cols()+2*C.num_ghost);
    extrapCopyCoords(xc_g, yc_g, xc, yc, C); // get ghost cell coord


    //i or x dir areas the columns of the matrix
    RowMajorMatrixXd Ai(xc.rows(), xn.cols()); // sets (row,col) (j,i)
    //j or y dir areas the rows or the matrix
    RowMajorMatrixXd Aj(xn.rows(), xc.cols()); // sets (row,col) (j,i)
    RowMajorMatrixXd Volume(xc.rows(), xc.cols());

    // SETUP NORMAL VECS
    //Ai normal vector has x and y components (so same size as Ai)
    RowMajorMatrixXd nix(Ai.rows(), Ai.cols());
    RowMajorMatrixXd niy(Ai.rows(), Ai.cols());
    //Aj normal vector has x and y componenets
    RowMajorMatrixXd njx(Aj.rows(), Aj.cols());
    RowMajorMatrixXd njy(Aj.rows(), Aj.cols());

    computeArea(Ai, Aj, xn, yn); // take nodal coords extract interface A // Checked
    computeNormalVectors(nix, niy, njx, njy,
        xn, yn, Ai, Aj); // grab the norm vec 4 computational domain // Checked
    computeVolume(Volume, xn, yn); // Checked

    outputArray(outputFolder, "xc_g", xc_g, 0);
    outputArray(outputFolder, "yc_g", yc_g, 0);
    outputArray(outputFolder, "xc", xc, 0);
    outputArray(outputFolder, "yc", yc, 0);

    outputArray(outputFolder, "Ai", Ai, 0);
    outputArray(outputFolder, "Aj", Aj, 0);
    outputArray(outputFolder, "njx", njx, 0);
    outputArray(outputFolder, "njy", njy, 0);
    outputArray(outputFolder, "nix", nix, 0);
    outputArray(outputFolder, "niy", niy, 0);
    outputArray(outputFolder, "volume", Volume, 0);

    // Now we can send the information that the different blocks need to 

    cout << " Prep to send " << endl;
    int nx = xc.cols();
    int ny = xc.rows();

    double *xn1d = xn.data();
    double *yn1d = yn.data();

    int nxb = nx/ni;
    int nyb = ny/nj;

    int nxbg = nxb+C.num_ghost*2;
    int nybg = nyb+C.num_ghost*2;

    /* This is coord for rank 0 or block 0
       printf("nxbg = %d, nybg = %d, nx = %d, ny = %d, nxb = %d, nyb = %d\n", nxbg, nybg, nx, ny, nxb, nyb);
       cout << xc_g.block(0, 0, nybg-1, nxbg-1) << endl;
       */

    nb[0]=nxb;
    nb[1]=nyb;
    for(int r = 1; r < size; r++)
      MPI_Send(nb, 2, MPI_INT, r, 999, com2d);

    /* All the variables geometry wise needed locally to loop
    RowMajorMatrixXd xcL_g, ycL_g;
    RowMajorMatrixXd nixL, niyL;
    RowMajorMatrixXd njxL, njyL;
    RowMajorMatrixXd AjL, AiL;
    RowMajorMatrixXd VolumeL;
    */

    // Variables resize to send
    xcL_g.resize(nybg, nxbg), ycL_g.resize(nybg, nxbg);
    nixL.resize(nyb, nxb+1), niyL.resize(nyb, nxb+1);
    njxL.resize(nyb+1, nxb), njyL.resize(nyb+1, nxb);
    AiL.resize(nyb, nxb+1), AjL.resize(nyb+1, nxb);
    VolumeL.resize(nyb, nxb);
    RowMajorMatrixXd temp;
    for(int j = 0; j < nj; j++)
      for(int i = 0; i < ni; i++)
      {
        int r = j*ni+i;
        if (r >= 1)
        {
          // send xcg, and ycg
          temp = xc_g.block(nyb*j, nxb*i, nybg, nxbg);
          MPI_Send(temp.data(), nxbg*nybg, MPI_DOUBLE, r, 111, com2d);
          temp = yc_g.block(nyb*j, nxb*i, nybg, nxbg);
          MPI_Send(temp.data(), nxbg*nybg, MPI_DOUBLE, r, 222, com2d);

          // send nix, and niy
          temp = nix.block(nyb*j, nxb*i, nyb, nxb+1);
          MPI_Send(temp.data(), nyb*(nxb+1), MPI_DOUBLE, r, 333, com2d);
          temp = niy.block(nyb*j, nxb*i, nyb, nxb+1);
          MPI_Send(temp.data(), nyb*(nxb+1), MPI_DOUBLE, r, 444, com2d);

          // send njx, and njy
          temp = njx.block(nyb*j, nxb*i, nyb+1, nxb);
          MPI_Send(temp.data(), nxb*(nyb+1), MPI_DOUBLE, r, 555, com2d);
          temp = njy.block(nyb*j, nxb*i, nyb+1, nxb);
          MPI_Send(temp.data(), nxb*(nyb+1), MPI_DOUBLE, r, 666, com2d);


          // send Aj, and Ai
          temp = Ai.block(nyb*j, nxb*i, nyb, nxb+1);
          MPI_Send(temp.data(), nyb*(nxb+1), MPI_DOUBLE, r, 777, com2d);
          temp = Aj.block(nyb*j, nxb*i, nyb+1, nxb);
          MPI_Send(temp.data(), nxb*(nyb+1), MPI_DOUBLE, r, 888, com2d);

          // send Volume
          temp = Volume.block(nyb*j, nxb*i, nyb, nxb);
          MPI_Send(temp.data(), nyb*nxb, MPI_DOUBLE, r, 1111, com2d);
        }
        else // set your own coords
        {
          xcL_g = xc_g.block(nyb*j, nxb*i, nybg, nxbg);
          ycL_g = yc_g.block(nyb*j, nxb*i, nybg, nxbg);

          nixL = nix.block(nyb*j, nxb*i, nyb, nxb+1);
          niyL = niy.block(nyb*j, nxb*i, nyb, nxb+1);

          njxL = njx.block(nyb*j, nxb*i, nyb+1, nxb);
          njyL = njy.block(nyb*j, nxb*i, nyb+1, nxb);

          AiL = Ai.block(nyb*j, nxb*i, nyb, nxb+1);
          AjL = Aj.block(nyb*j, nxb*i, nyb+1, nxb);

          VolumeL = Volume.block(nyb*j, nxb*i, nyb, nxb);
        }
      }
  }
  else
  {
    int source = 0;
    MPI_Status status;
    // receive sizes
    MPI_Recv(&nb, 2, MPI_INT, source, 999, com2d, &status);

    int nxb = nb[0];
    int nyb = nb[1];

    int nxbg = nxb+2*C.num_ghost;
    int nybg = nyb+2*C.num_ghost;

    xcL_g.resize(nybg, nxbg);
    ycL_g.resize(nybg, nxbg);
    double *xcL_gP = new double [(nxbg)*(nybg)];
    double *ycL_gP = new double [(nxbg)*(nybg)];
    // get xcg and ycg locall for ranks 1->size
    MPI_Recv(xcL_gP, nxbg*nybg, MPI_DOUBLE, source, 111, com2d, &status);
    MPI_Recv(ycL_gP, nxbg*nybg, MPI_DOUBLE, source, 222, com2d, &status);
    for(int j = 0; j < nybg; j++)
      for(int i = 0; i < nxbg; i++)
      {
        xcL_g(j,i) = xcL_gP[j*nxbg+i];
        ycL_g(j,i) = ycL_gP[j*nxbg+i];
      }
    delete[] xcL_gP; xcL_gP=NULL;
    delete[] ycL_gP; ycL_gP=NULL;

    nixL.resize(nyb, nxb+1), niyL.resize(nyb, nxb+1);
    double *nixLP = new double [nyb*(nxb+1)];
    double *niyLP = new double [nyb*(nxb+1)];
    MPI_Recv(nixLP, nyb*(nxb+1), MPI_DOUBLE, source, 333, com2d, &status);
    MPI_Recv(niyLP, nyb*(nxb+1), MPI_DOUBLE, source, 444, com2d, &status);
    for(int j = 0; j < nyb; j++)
      for(int i = 0; i < nxb+1; i++)
      {
        nixL(j,i) = nixLP[j*(nxb+1)+i];
        niyL(j,i) = niyLP[j*(nxb+1)+i];
      }
    delete[] nixLP; nixLP=NULL;
    delete[] niyLP; niyLP=NULL;


    njxL.resize(nyb+1, nxb), njyL.resize(nyb+1, nxb);
    double *njxLP = new double [nxb*(nyb+1)];
    double *njyLP = new double [nxb*(nyb+1)];
    MPI_Recv(njxLP, nxb*(nyb+1), MPI_DOUBLE, source, 555, com2d, &status);
    MPI_Recv(njyLP, nxb*(nyb+1), MPI_DOUBLE, source, 666, com2d, &status);
    for(int j = 0; j < nyb+1; j++)
      for(int i = 0; i < nxb; i++)
      {
        njxL(j,i) = njxLP[j*nxb+i];
        njyL(j,i) = njyLP[j*nxb+i];
      }
    delete[] njxLP; njxLP=NULL;
    delete[] njyLP; njyLP=NULL;


    AiL.resize(nyb, nxb+1), AjL.resize(nyb+1, nxb);
    double *AiLP = new double [nyb*(nxb+1)];
    double *AjLP = new double [nxb*(nyb+1)];

    MPI_Recv(AiLP, nyb*(nxb+1), MPI_DOUBLE, source, 777, com2d, &status);
    MPI_Recv(AjLP, nxb*(nyb+1), MPI_DOUBLE, source, 888, com2d, &status);
    for(int j = 0; j < nyb; j++)
      for(int i = 0; i < nxb+1; i++)
        AiL(j,i) = AiLP[j*(nxb+1)+i];
    for(int j = 0; j < nyb+1; j++)
     for(int i = 0; i < nxb; i++)
        AjL(j,i) = AjLP[j*(nxb)+i];
    delete[] AiLP; AiLP=NULL;
    delete[] AjLP; AjLP=NULL;


    VolumeL.resize(nyb, nxb);
    double *VolumeLP = new double [nyb*nxb];
    MPI_Recv(VolumeLP, nyb*nxb, MPI_DOUBLE, source, 1111, com2d, &status);
    for(int j = 0; j < nyb; j++)
     for(int i = 0; i < nxb; i++)
        VolumeL(j,i) = VolumeLP[j*nxb+i];
 
  }
    
  outputArray(outputFolder, "xcL_g", xcL_g, rank);
  outputArray(outputFolder, "ycL_g", ycL_g, rank);
  outputArray(outputFolder, "nixL", nixL, rank);
  outputArray(outputFolder, "niyL", niyL, rank);
  outputArray(outputFolder, "njxL", njxL, rank);
  outputArray(outputFolder, "njyL", njyL, rank);
  outputArray(outputFolder, "AiL", AiL, rank);
  outputArray(outputFolder, "AjL", AjL, rank);
 

  return com2d;

  //printf("\nmy rank = %d, and I have nxc = %d, and nyb = %d\n", rank, nb[0], nb[1]);
}

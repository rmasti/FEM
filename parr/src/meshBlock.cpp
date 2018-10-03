/*
 * meshBlock.cpp
 * Written By: Robert Masti
 * 10/03/2018
 * This function has rank 0 compute areas, volume, coord, normal vecs
 * and distributes them to the other ranks
 * It also returns the cart 2d communicator
 */
#include "mhdRT.hpp"

MPI_Comm meshBlock(
    const string mesh,            //In: mesh name to read
    const string outputFolder,    //In: name of output dir
    RowMajorMatrixXd &xcL_g,      //Out: local x center coord
    RowMajorMatrixXd &ycL_g,      //Out: local y center coord
    RowMajorMatrixXd &nixL,       //Out: local x comp norm vec i dir
    RowMajorMatrixXd &niyL,       //Out: local y comp norm vec i dir
    RowMajorMatrixXd &njxL,       //Out: local x comp norm vec j dir
    RowMajorMatrixXd &njyL,       //Out: local y comp norm vec j dir
    RowMajorMatrixXd &AiL,        //Out: local area i dir
    RowMajorMatrixXd &AjL,        //Out: local area j dir
    RowMajorMatrixXd &VolumeL,    //Out: local cell volume
    const constants C)                  
{
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size); // number of proc
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // current proc

  MPI_Comm com2d; // Setup cartesian communicator

  int ni, nj; 
  int ng = C.num_ghost;

  int cnt=size; // cnt will give n of 2^n
  int n=0;
  while(cnt >=2)
  {
    cnt = cnt/2;
    n++;
  }
  nj = pow(2, (n/2));         // proc vert dir
  ni = pow(2, (n/2) + n%2);   // proc horz dir

  int dim[2] = {nj, ni};      // 2 int array for MPI_Cart_Create
  int reorder = FALSE;        // ReOrder for fastest comm
  int period[2];              // periodic boundary partit
  period[1] = TRUE;            
  period[0] = FALSE;
  int coord[2];               // current proc coord
  int id;                      

  int u, d, l, r; // up down left right

  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &com2d);
  MPI_Comm_rank(com2d, &rank);
  int nb[2];

  // have the first processor read in data and give it to all other ranks
  if(rank == 0)
  {
    RowMajorMatrixXd xn;
    RowMajorMatrixXd yn;
    RowMajorMatrixXd xc, yc;
    // Read in mesh 
    inputMesh(xn,yn,xc,yc, mesh);

    // extrap to ghost cells
    RowMajorMatrixXd xc_g(xc.rows()+2*ng, xc.cols()+2*ng);
    RowMajorMatrixXd yc_g(yc.rows()+2*ng, yc.cols()+2*ng);
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

    // fill matrices
    computeArea(Ai, Aj, xn, yn); 
    computeNormalVectors(nix, niy, njx, njy, xn, yn, Ai, Aj); 
    computeVolume(Volume, xn, yn);

    // output for plotting purposes
    outputArray(outputFolder, "xc_g", xc_g, 0);
    outputArray(outputFolder, "yc_g", yc_g, 0);
    outputArray(outputFolder, "xc", xc, 0);
    outputArray(outputFolder, "yc", yc, 0);
    outputArray(outputFolder, "xn", xn, 0);
    outputArray(outputFolder, "yn", yn, 0);
    outputArray(outputFolder, "Ai", Ai, 0);
    outputArray(outputFolder, "Aj", Aj, 0);
    outputArray(outputFolder, "njx", njx, 0);
    outputArray(outputFolder, "njy", njy, 0);
    outputArray(outputFolder, "nix", nix, 0);
    outputArray(outputFolder, "niy", niy, 0);
    outputArray(outputFolder, "volume", Volume, 0);

    // Send the information of sizes to other ranks
    int nx = xc.cols();
    int ny = xc.rows();

    int nxb = nx/ni; // nx per block
    int nyb = ny/nj;

    int nxbg = nxb+ng*2;
    int nybg = nyb+ng*2;

    nb[0]=nxb;
    nb[1]=nyb;
    for(int r = 1; r < size; r++)
      MPI_Send(nb, 2, MPI_INT, r, 999, com2d);

    // Variables resize to send
    xcL_g.resize(nybg, nxbg), ycL_g.resize(nybg, nxbg);
    nixL.resize(nyb, nxb+1), niyL.resize(nyb, nxb+1);
    njxL.resize(nyb+1, nxb), njyL.resize(nyb+1, nxb);
    AiL.resize(nyb, nxb+1), AjL.resize(nyb+1, nxb);
    VolumeL.resize(nyb, nxb);
    RowMajorMatrixXd temp;
    for(int j = 0; j < nj; j++)
    {
      for(int i = 0; i < ni; i++)
      {
        int r = j*ni+i;
        if (r >= 1) // send data to other ranks else
        {
          // send xcg, and ycg
          temp = xc_g.block(nyb*j, nxb*i, nybg, nxbg); // grab block 
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
  }
  else
  {
    int source = 0; // rank 0 sent all 
    MPI_Status status;
    // receive sizes
    MPI_Recv(&nb, 2, MPI_INT, source, 999, com2d, &status);

    int nxb = nb[0];
    int nyb = nb[1];

    int nxbg = nxb+2*ng;
    int nybg = nyb+2*ng;


    // get xcg and ycg locall for ranks 1->size
    xcL_g.resize(nybg, nxbg), ycL_g.resize(nybg, nxbg);
    MPI_Recv(xcL_g.data(), nxbg*nybg, MPI_DOUBLE, source, 111, com2d, &status);
    MPI_Recv(ycL_g.data(), nxbg*nybg, MPI_DOUBLE, source, 222, com2d, &status);
    
    // fill normal vecs i dir
    nixL.resize(nyb, nxb+1), niyL.resize(nyb, nxb+1);
    MPI_Recv(nixL.data(), nyb*(nxb+1), MPI_DOUBLE, source, 333, com2d, &status);
    MPI_Recv(niyL.data(), nyb*(nxb+1), MPI_DOUBLE, source, 444, com2d, &status);

    // fill normal vecs j dir
    njxL.resize(nyb+1, nxb), njyL.resize(nyb+1, nxb);
    MPI_Recv(njxL.data(), nxb*(nyb+1), MPI_DOUBLE, source, 555, com2d, &status);
    MPI_Recv(njyL.data(), nxb*(nyb+1), MPI_DOUBLE, source, 666, com2d, &status);

    // fill areas
    AiL.resize(nyb, nxb+1), AjL.resize(nyb+1, nxb);
    MPI_Recv(AiL.data(), nyb*(nxb+1), MPI_DOUBLE, source, 777, com2d, &status);
    MPI_Recv(AjL.data(), nxb*(nyb+1), MPI_DOUBLE, source, 888, com2d, &status);

    // fill volumes
    VolumeL.resize(nyb, nxb);
    MPI_Recv(VolumeL.data(), nyb*nxb, MPI_DOUBLE, source, 1111, com2d, &status);
  }

  outputArray(outputFolder, "xcL_g", xcL_g, rank);
  outputArray(outputFolder, "ycL_g", ycL_g, rank);
  outputArray(outputFolder, "nixL", nixL, rank);
  outputArray(outputFolder, "niyL", niyL, rank);
  outputArray(outputFolder, "njxL", njxL, rank);
  outputArray(outputFolder, "njyL", njyL, rank);
  outputArray(outputFolder, "AiL", AiL, rank);
  outputArray(outputFolder, "AjL", AjL, rank);


  MPI_Barrier(com2d);
  return com2d;
}

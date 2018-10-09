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
    RowMajorMatrixXd &xcLg,      //Out: local x center coord
    RowMajorMatrixXd &ycLg,      //Out: local y center coord
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

  int dim[2] = {nj, ni};        // 2 int array for MPI_Cart_Create
  int reorder = FALSE;          // Reorder for fastest comm
  int period[2] = {FALSE, TRUE};// periodic boundary partit

  MPI_Comm com2d; // Setup cartesian communicator
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &com2d);
  MPI_Comm_rank(com2d, &rank);

  int nb[2];
  // have the first processor read in data and give it to all other ranks
  if(rank == 0)
  {
    RowMajorMatrixXd xn, yn;
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
    xcLg.resize(nybg, nxbg), ycLg.resize(nybg, nxbg);
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
          xcLg = xc_g.block(nyb*j, nxb*i, nybg, nxbg);
          ycLg = yc_g.block(nyb*j, nxb*i, nybg, nxbg);

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
    xcLg.resize(nybg, nxbg), ycLg.resize(nybg, nxbg);
    MPI_Recv(xcLg.data(), nxbg*nybg, MPI_DOUBLE, source, 111, com2d, &status);
    MPI_Recv(ycLg.data(), nxbg*nybg, MPI_DOUBLE, source, 222, com2d, &status);
    
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

  outputArray(outputFolder, "xcLg", xcLg, rank);
  outputArray(outputFolder, "ycLg", ycLg, rank);
  outputArray(outputFolder, "nixL", nixL, rank);
  outputArray(outputFolder, "niyL", niyL, rank);
  outputArray(outputFolder, "njxL", njxL, rank);
  outputArray(outputFolder, "njyL", njyL, rank);
  outputArray(outputFolder, "AiL", AiL, rank);
  outputArray(outputFolder, "AjL", AjL, rank);


  MPI_Barrier(com2d);
  return com2d;
}

void inputMesh(
    RowMajorMatrixXd& xn,         // Output - x coord nodes
    RowMajorMatrixXd& yn,         // Output - y coord nodes
    RowMajorMatrixXd& xc,         // Output - x coord centers
    RowMajorMatrixXd& yc,         // Output - y coord centers
    const string mesh)    // input - mesh location
{
  ifstream infile;
  infile.open(mesh.c_str());

  if (!infile){
    cerr << "ERROR: Unable to open the mesh file" << endl;
    exit(1);
  }

  double readval;
  string line;
  int vals[3];
  for(int i = 0; i < 3; i++ )
  {
    getline(infile, line);
    istringstream streamA(line);
    streamA >> readval;
    vals[i] = readval;
  }
  int ntot = vals[0];
  int nx = vals[1];
  int ny = vals[2];

  //double dat[ntot*2];
  VectorXd dat(ntot*2);
  int i = 0;
  while (infile.good())
    while (getline(infile, line))
    {
      istringstream streamA(line);
      while(streamA >> readval)
        dat(i) = readval;
      i++;
    } 

  xn.resize(ny, nx); //number of rows, number of columns
  yn.resize(ny, nx); //number of rows, number of columns
  xc.resize(ny-1, nx-1);
  yc.resize(ny-1, nx-1);

  int offset = ntot*2;
  for(int j = 0; j < ny; j++)
    for(int i = 0; i < nx; i++)
    {
      int xInd = j*nx+i;
      int yInd = j*nx+i+nx*ny;
      xn(j, i) = dat[xInd];
      yn(j, i) = dat[yInd];
    }
  for (int j = 0; j < xc.rows(); j++)
    for (int i = 0 ; i < xc.cols(); i++)
    {
      // take average of all four xvals from corn
      // order: botL + topL + botR + topR 
      xc(j,i) = 0.25*(xn(j,i) + xn(j+1,i) + 
          xn(j,i+1) + xn(j+1,i+1));
      yc(j,i) = 0.25*(yn(j,i) + yn(j+1,i) + 
          yn(j,i+1) + yn(j+1,i+1));
    }
}

void extrapCopyCoords(
    // This function will copy the interior coord of the cells, but
    // will also extrapolate those coords to the ghost cells
    RowMajorMatrixXd& xc_g,        // output - cell x coord with ghosts 
    RowMajorMatrixXd& yc_g,        // output - cell y coord with ghosts
    const RowMajorMatrixXd& xc,    // input - cell x coord without ghosts
    const RowMajorMatrixXd& yc,    // input - cell y coord without ghosts
    constants C            // input - constants C for num ghost
    )
{
  int ni = xc.cols();
  int nj = xc.rows();
  // copy the interior values for both x and y
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      int i_g = i+C.num_ghost; // set iter to 3+i
      int j_g = j+C.num_ghost; // set iter to 3+j
      xc_g(j_g, i_g) = xc(j,i);
      yc_g(j_g, i_g) = yc(j,i);
    }
  }
  int bot = C.num_ghost - 1; // 1st ghost cell based off 
  int top = C.num_ghost + nj; // top index start
  int left = C.num_ghost -1; // same as bot
  int right = C.num_ghost + ni; // 1st ghost on right
  // extrapolate to the ghost cells now

  // j dir which means the interior need extrap
  for (int i = C.num_ghost; i < ni+C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bottom index - j is the filling of ghost
      xc_g(bot-j,i) = 2.0*xc_g(bot-j+1,i) - xc_g(bot-j+2,i);
      yc_g(bot-j,i) = 2.0*yc_g(bot-j+1,i) - yc_g(bot-j+2,i);
      // top index top + j is filling ghost
      xc_g(top+j,i) = 2.0*xc_g(top+j-1,i) - xc_g(top+j-2,i);
      yc_g(top+j,i) = 2.0*yc_g(top+j-1,i) - yc_g(top+j-2,i);
    }
  }
  // i dir which means the interior need extrap
  for (int i = 0; i < C.num_ghost; i ++)
  {
    for (int j = C.num_ghost; j < C.num_ghost+nj; j++)
    {
      // left index - i is the filling of ghost
      xc_g(j,left-i) = 2.0*xc_g(j,left-i+1) - xc_g(j,left-i+2);
      yc_g(j,left-i) = 2.0*yc_g(j,left-i+1) - yc_g(j,left-i+2);
      // right index + i is the filling of ghost
      xc_g(j,right+i) = 2.0*xc_g(j,right+i-1) - xc_g(j,right+i-2);
      yc_g(j,right+i) = 2.0*yc_g(j,right+i-1) - yc_g(j,right+i-2);
    }
  }
  //set to 1 out the corners for sake of clarity.
  for (int i = 0; i < C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bot left
      xc_g(j,i) = 1.0; yc_g(j,i) = 1.0;
      // bot right
      xc_g(j,right+i) = 1.0; yc_g(j,right+i) = 1.0;
      // top right
      xc_g(top+j,right+i) = 1.0; yc_g(top+j,right+i) = 1.0;
      // top left
      xc_g(top+j,i) = 1.0; yc_g(top+j,i) = 1.0;
    }
  }
}

void computeArea(
    // This function computes the area based on the node locations 
    // for both the i facing direction and the j facing direction
    RowMajorMatrixXd& Ai,          // output - i-direction areas (c, n)size    
    RowMajorMatrixXd& Aj,          // output - j-direction areas (n, c)size
    const RowMajorMatrixXd& xn,    // input - x coord for nodes 
    const RowMajorMatrixXd& yn     // input - y coord for nodes
    )
{
  // loop for Ai
  // NOTE i is the columns and j is the rows
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      double dx = xn(j,i) - xn(j+1,i); // compute bottom - top x
      double dy = yn(j,i) - yn(j+1,i); // compute bottom - top y
      Ai(j,i) = sqrt( dx*dx + dy*dy ); // pythag
    }
  }
  // loop for Aj
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      double dx = xn(j,i) - xn(j,i+1); // compute left - right x
      double dy = yn(j,i) - yn(j,i+1); // compute left - right y
      Aj(j,i) = sqrt( dx*dx + dy*dy);  // pytha
    }
  }
  // NOTE area is the dist times a w into the z dir which is just 1 for 2d
}

void computeNormalVectors(
    // This function takes in the previously calculated areas and
    // Finds the fraction in the x and y physical directions
    RowMajorMatrixXd& nix,    // output - i with A of normal vector in x phys
    RowMajorMatrixXd& niy,    // output - i comp but the yhat phys
    RowMajorMatrixXd& njx,    // output - j with A of normal vector in x phys
    RowMajorMatrixXd& njy,    // output - j dir but with yhat phys
    const RowMajorMatrixXd& xn,    // input - nodal x coordinate
    const RowMajorMatrixXd& yn,    // input - nodal y coordinate
    const RowMajorMatrixXd& Ai,    // input - Area pointed in i (x square grid)
    const RowMajorMatrixXd& Aj     // input - Area pointed in j (y square grid)
    )
{
  // NOTE i is the columns and j are the rows for A(row,col)
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      // See notes Sec 6 slide 75
      nix(j,i) = (yn(j+1,i) - yn(j,i)) / Ai(j,i); 
      niy(j,i) = -(xn(j+1,i) - xn(j,i)) / Ai(j,i);
    }
  }
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      njx(j,i) = -(yn(j,i+1) - yn(j,i)) / Aj(j,i);
      njy(j,i) = (xn(j,i+1) - xn(j,i)) / Aj(j,i);
    }
  }
}

void computeVolume(
    // This function computes the volumes of the interior cell points. 
    RowMajorMatrixXd& Volume,      // output - Volume of the cells
    const RowMajorMatrixXd& xn,    // input - x coord for nodal vals
    const RowMajorMatrixXd& yn     // input - y coord for nodal vals
    )
{
  double dx1, dx2, dy1, dy2; // find the diagonals on the cell
  double cross; // cross product
  for (int i = 0; i < Volume.cols(); i++)
  {
    for (int j = 0; j < Volume.rows(); j++)
    {
      dx1 = xn(j,i) - xn(j+1,i+1); // bot L - top R
      dy1 = yn(j,i) - yn(j+1,i+1); // bot L - top R
      dx2 = xn(j+1,i) - xn(j,i+1); // top L - bot R
      dy2 = yn(j+1,i) - yn(j,i+1); // top L - bot R
      cross = abs(dx1*dy2 - dx2*dy1); // cross product
      Volume(j,i) = 0.5*cross; // assume w=1 volume of cell
    }
  } 
}

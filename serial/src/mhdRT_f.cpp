# include "mhdRT.hpp"


void extrapCopyCoords(
    double xl_g[], 
    double xr_g[], 
    double xb_g[], 
    double xt_g[],  
    double yl_g[], 
    double yr_g[], 
    double yb_g[], 
    double yt_g[],  
    double xc[], 
    double yc[], 
    constants C)
{
  if(C.num_ghost < 2)
  {
    cerr << "ERROR: Need to use 2 or more Ghost Layers" << endl;
    exit(-1);
  }
  for( int i=0; i < C.ny_c; i++)
  {
    // Left 
    xl_g[i*C.num_ghost] = 2.0*xc[i*C.nx_c] - xc[i*C.nx_c+1];
    xl_g[i*C.num_ghost+1] = 2.0*xl_g[i*C.num_ghost] - xc[i*C.nx_c];
    yl_g[i*C.num_ghost] = 2.0*yc[i*C.nx_c] - yc[i*C.nx_c+1];
    yl_g[i*C.num_ghost+1] = 2.0*yl_g[i*C.num_ghost] - yc[i*C.nx_c];

    // Right
    xr_g[i*C.num_ghost] = 2.0*xc[(i+1)*(C.nx_c)-1] - xc[(i+1)*(C.nx_c)-2];
    xr_g[i*C.num_ghost+1] = 2.0*xr_g[i*C.num_ghost] - xc[(i+1)*(C.nx_c)-1];
    yr_g[i*C.num_ghost] = 2.0*yc[(i+1)*(C.nx_c)-1] - yc[(i+1)*(C.nx_c)-2];
    yr_g[i*C.num_ghost+1] = 2.0*yr_g[i*C.num_ghost] - yc[(i+1)*(C.nx_c)-1];

//printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xr_g[i*C.num_ghost+1],xr_g[i*C.num_ghost] , xc[(i+1)*(C.nx_c)-1], xc[(i+1)*(C.nx_c)-2] );

    // Extrapolate the rest 
    for( int j=2; j < C.num_ghost; j++)
    {
      xr_g[i*C.num_ghost+j] = 2.0*xr_g[i*C.num_ghost+j-1] -xr_g[i*C.num_ghost+j-2];
      yr_g[i*C.num_ghost+j] = 2.0*yr_g[i*C.num_ghost+j-1] -yr_g[i*C.num_ghost+j-2];
      xl_g[i*C.num_ghost+j] = 2.0*xl_g[i*C.num_ghost+j-1] -xl_g[i*C.num_ghost+j-2];
      yl_g[i*C.num_ghost+j] = 2.0*yl_g[i*C.num_ghost+j-1] -yl_g[i*C.num_ghost+j-2];
    }
  }
  // top bottom
  for(int j=0; j < C.nx_c; j++)
  {
    // top 
    xt_g[j*C.num_ghost] = 2.0*xc[C.nx_c*(C.ny_c-1)+j] - xc[C.nx_c*(C.ny_c-2)+j];
    xt_g[j*C.num_ghost+1] = 2.0*xt_g[j*C.num_ghost] - xc[C.nx_c*(C.ny_c-1)+j] ;
    yt_g[j*C.num_ghost] = 2.0*yc[C.nx_c*(C.ny_c-1)+j] - yc[C.nx_c*(C.ny_c-2)+j];
    yt_g[j*C.num_ghost+1] = 2.0*yt_g[j*C.num_ghost] - yc[C.nx_c*(C.ny_c-1)+j] ;

    // bottom 
    xb_g[j*C.num_ghost] = 2.0*xc[j] - xc[C.nx_c+j];
    xb_g[j*C.num_ghost+1] = 2.0*xb_g[j*C.num_ghost] - xc[j] ;
    yb_g[j*C.num_ghost] = 2.0*yc[j] - yc[C.nx_c+j];
    yb_g[j*C.num_ghost+1] = 2.0*yb_g[j*C.num_ghost] - yc[j] ;
 
    //printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",yb_g[j*C.num_ghost+1],yb_g[j*C.num_ghost] , yc[j], yc[C.nx_c+j] );
    //printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xb_g[j*C.num_ghost+1],xb_g[j*C.num_ghost] , xc[C.nx_c*(C.ny_c-1)+j], xc[C.nx_c*(C.ny_c-2)+j] );
    //printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n", j*C.num_ghost,C.nx_c*(C.ny_c-1)+j ,C.nx_c*(C.ny_c-2)+j) ;

    for(int i=2; i < C.num_ghost; i++)
    {
      xt_g[j*C.num_ghost+i] = 2.0*xt_g[j*C.num_ghost+i-1] -xt_g[j*C.num_ghost+i-2];
      yt_g[j*C.num_ghost+i] = 2.0*yt_g[j*C.num_ghost+i-1] -yt_g[j*C.num_ghost+i-2];
      xb_g[j*C.num_ghost+i] = 2.0*xb_g[j*C.num_ghost+i-1] -xb_g[j*C.num_ghost+i-2];
      yb_g[j*C.num_ghost+i] = 2.0*yb_g[j*C.num_ghost+i-1] -yb_g[j*C.num_ghost+i-2];
    }
  }
}


void getCoord(double xn[], double yn[], double xc[], double yc[], constants C)
{
  int nx = C.nx_c+1;
  int ny = C.ny_c+1;

  string filename = readMeshName(C);
  const char *FILENAME = filename.c_str();
  //open the instream
  ifstream infile;
  infile.open(FILENAME);

  double readval;
  vector<double> data;
  string line;

  if (infile.is_open())
  {
    while ( getline (infile,line))
    {
      istringstream streamA(line);
      while (streamA >> readval)
      {
        data.push_back(readval);
      }
    }
    infile.close();
  }

  for (int i=0; i < nx*ny; i++)
  {
    int datindex = i+3;
    xn[i] = data[datindex];
    yn[i] = data[datindex + nx*ny];
  }

  for (int i=0; i < C.ny_c; i++)
  {
    for (int j=0; j < C.nx_c; j++)
    {

      xc[i*C.nx_c+j] = 0.25*( xn[nx*i+j] + xn[nx*i+j+1] + xn[nx*(i+1)+j] + xn[nx*(i+1)+j+1] );
      yc[i*C.nx_c+j] = 0.25*( yn[nx*i+j] + yn[nx*i+j+1] + yn[nx*(i+1)+j] + yn[nx*(i+1)+j+1] );

      //cout << nx*i+j <<" " << nx*i+j+1  <<  " " << nx*(i+1)+j  << " " <<nx*(i+1)+j+1  << endl;
      //cout << xc[i*C.nx_c+j] <<" " << yc[i*C.nx_c+j] <<  " " << xn[nx*i+j] << endl;
    }
  }
}

void meshSize(int* nx_i, int* ny_i, constants C)
{
  string filename = readMeshName(C);
  const char *FILENAME = filename.c_str();
  //open the instream
  ifstream infile;
  infile.open(FILENAME);

  int readval;
  vector<int> data;
  string line;
  if (infile.is_open())
  {
    int i=0;
    while ( getline (infile,line) && i < 3)
    {
      istringstream streamA(line);
      while (streamA >> readval)
      {
        data.push_back(readval);
      }
      i++;
    }
    infile.close();
  }

  *nx_i = data[1];
  *ny_i = data[2];

}

string readMeshName(
    // This function reads in the mesh name and return the filename
    // it will be used to extract the data
    constants C
    ) 
{   

  if (C.f_mesh == 1) // read in smallest mesh
    return "mesh/debugMatlab.msh";
  if (C.f_mesh == 2) 
    return "grids/Curvilinear/2d17.grd";

  return "failure"; // return failure because it can not recognize the file
}

    /*
       cout << C.nx_c << " " << C.ny_c<< endl;
       printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xr_g[i*C.num_ghost+1],xr_g[i*C.num_ghost] , xc[(i+1)*(C.nx_c)-1], xc[(i+1)*(C.nx_c)-2] );
       printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n",i*C.num_ghost ,  (i+1)*C.nx_c-1, (i+1)*C.nx_c-2);
       printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n",i*C.num_ghost ,  i*C.nx_c, i*C.nx_c+1);
       printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xl_g[i*C.num_ghost+1],xl_g[i*C.num_ghost] ,  xc[i*C.nx_c], xc[i*C.nx_c+1]);
       printf("xc1stInt = %lf, and xc2ndInt = %lf \n", xc[i*C.ny_c], xc[i*C.ny_c+1]);
       xl_g[i*C.num_ghost] = 2.0*xc[i*C.nx_c+1] - xc[i*C.nx_c+2];
       xl_g[i*C.num_ghost+1] = 2.0*xl_g[i*C.num_ghost] - xc[i*C.nx_c+1];

       printf("xl1st = %lf, and xl2nd = %lf \n", xl_g[i*C.num_ghost], xl_g[i*C.num_ghost+1]);
    */


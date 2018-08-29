# include "mhdRT.hpp"


void cons2Prim(double V[], double U[], int size)
{
  int index; 
  double U0[NEQ];
  for(int i = 0; i < size/NEQ; i++)
  {
    index = i*NEQ;
    for(int eq = 0; eq < NEQ; eq++)
      U0[eq] = U[index+eq], V[index+eq] = U0[eq];

    V[index+uid] = U0[uid]/U0[rhoid]; //rhou/rho
    V[index+vid] = U0[vid]/U0[rhoid]; //rhou/rho
    V[index+wid] = U0[wid]/U0[rhoid]; //rhou/rho
    V[index+pid] = (GAMMA-1.0)*(U0[pid] - 0.5*(U0[uid]*U0[uid] + U0[vid]*U0[vid] + U0[wid]*U0[wid])/U0[rhoid] - 0.5*((U0[bxid]*U0[bxid] + U0[byid]*U0[byid] + U0[bzid]*U0[bzid]))/MU0);
  }
}

void prim2Cons(double U[], double V[], int size)
{
  int index;
  double V0[NEQ];
  for(int i = 0; i < size/NEQ; i++)
  {
    index = i*NEQ;
    for(int eq = 0; eq < NEQ; eq++)
      V0[eq] = V[index+eq], U[index+eq] = V0[eq];
  
    U[index+uid] = V0[rhoid]*V0[uid]; // rhou
    U[index+vid] = V0[rhoid]*V0[vid]; // rhov
    U[index+wid] = V0[rhoid]*V0[wid]; // rhow
    U[index+pid] = V0[pid]/(GAMMA-1.0) + 0.5*V0[rhoid]*(V0[uid]*V0[uid] + V0[vid]*V0[vid] + V0[wid]*V0[wid]) + 0.5*((V0[bxid]*V0[bxid] + V0[byid]*V0[byid] + V0[bzid]*V0[bzid]))/MU0; // e 
  }
}

void initialize(double V[], double rc[], double thetc[], constants C)
{
  double V0[NEQ];
  double r, thet;
  int cIndex;
  switch (C.f_case)
  {
    case 1:
      double rhoL = 1.0;
      double rhoH = 2.0;
      double p0 = 2.5;
      double g = 0.1;
      double vPert = 0.01;
      double bparr0 = 0.0*sqrt(4.0*PI);
      double circ = 2.0*PI*0.375;
      double lambda = circ/20;
      double randPert;
      V0[wid] = 0.0;

      for(int i = 0; i < C.ny_c; i++)
      {
        for(int j = 0; j < C.nx_c; j++)
        {
          cIndex = i*C.nx_c + j;
          r = rc[cIndex];
          thet = thetc[cIndex];
          if(r < 0.375)
          {
            V0[rhoid] = rhoL;
            V0[pid] = p0-V0[rhoid]*g*(r-0.375);
          }
          else
          {
            V0[rhoid] = rhoH;
            V0[pid] = p0 - V[rhoid]*g*(r-0.375);
          }; 
          randPert = ((double) rand() / (RAND_MAX));
          double decay =exp(-20.0*(r-0.375)*(r-0.375)/(0.125*0.125)); 
          double thetVar = randPert*(2.0*PI/lambda)*thet+randPert*PI;
          V0[uid] = vPert*cos(thetVar)*decay; 
          V0[vid] = vPert*sin(thetVar)*decay; 
          V0[bxid] = bparr0*cos(thet);
          V0[byid] = bparr0*sin(thet);
          V0[bzid] = 0.0;
          /*
          for(int eq = 0; eq < NEQ; eq++)
            cout << V0[eq] << " "; 
          cout << endl;
          */
          // Fill in V
          for(int eq = 0; eq < NEQ; eq++)
            V[cIndex*NEQ+eq] = V0[eq];
            //V[i*C.nx_c+j*NEQ + eq]
        }
      }
  }
}

void computeVolume(double volume[], double xn[], double yn[], constants C)
{
  double dx1, dx2, dy1, dy2, cross;
  int nx = C.nx_c + 1;
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int cIndex = i*C.nx_c+j;
      int nIndex = i*nx+j;

      //bot L - top R
      dx1 = xn[nIndex] - xn[nIndex+nx+1];
      dy1 = yn[nIndex] - yn[nIndex+nx+1];
      //top L - bot R
      dx2 = xn[nIndex+nx] - xn[nIndex+1];
      dy2 = yn[nIndex+nx] - yn[nIndex+1];

      cross = abs(dx1*dy2 - dx2*dy1);
      volume[cIndex] = 0.5*cross;
    }
  }
}

void computeAreasAndNormalVectors(
    double njx[], 
    double njy[], 
    double nix[], 
    double niy[], 
    double Aj[], 
    double Ai[], 
    double xn[],
    double yn[], 
    constants C)
{
  int nx = C.nx_c+1; // number of interfaces/nodes x dir
  int ny = C.ny_c+1; // and y dir

  // loop for Aj or x dir
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      int nIndex = i*nx+j;
      double dx = xn[nIndex] - xn[nIndex+nx];
      double dy = yn[nIndex] - yn[nIndex+nx];
      Aj[nIndex] = sqrt(dx*dx + dy*dy);
      njx[nIndex] = dy / Aj[nIndex];
      njy[nIndex] = -1.0*dx / Aj[nIndex];
      //cout << nIndex << " " << nIndex+nx << endl;
    }
  }
  // loop for Ai or the ydir
  for(int i = 0; i < ny; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int nIndex = i*nx+j;
      double dx = xn[nIndex] - xn[nIndex+1];
      double dy = yn[nIndex+j] - yn[nIndex+1];

      int fIndex = i*C.nx_c+j;
      Ai[fIndex] = sqrt(dx*dx + dy*dy);    
      nix[fIndex] = -1.0*dy / Ai[fIndex];
      niy[fIndex] = dx / Ai[fIndex];
    }
  }
}

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


void computeRTheta(double rc[], double thetc[], double xc[], double yc[], constants C)
{
  int index; double x, y;
  for(int i = 0; i < C.ny_c; i++){
    for(int j = 0; j < C.nx_c; j++)
    {
      index = i*C.nx_c+j;
      x = xc[index]; y = yc[index];
      rc[index] = sqrt(x*x+y*y);
      thetc[index] = atan(y/x);
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

void outputArray(string Address, string FileName, double out[], int size,  int n)
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Suffix = ".txt";
  string Sub = "-";
  Address = Address + "/" + FileName + Sub + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component

  for(int i = 0; i < size; i++)
    outfile << setprecision(14) << out[i]  << endl;

  outfile.close();
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

//cout << Aj[i*nx+j] << endl;
//cout << dx <<" "  << dy  <<  endl;
//printf("faceIndex = %d, BotNode = %d,  TopNode = %d,\n", i*nx+j ,i*nx+j , (i+1)*nx+j) ;

//cout << dx << " " << dy << endl;
//cout << Ai[i*C.nx_c+j] << endl;
//printf("faceIndex = %d, leftNode = %d,  rightNode = %d,\n", i*C.nx_c+j, i*nx+j, i*nx+j+1 ) ;

*/


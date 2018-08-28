//# include "gmsh_io.hpp"
# include "mhdRT.hpp"


void extrapCopyCoords(double xc_g[], double yc_g[], double xc[], double yc[], constants C)
{

  int nx = C.nx_c;
  int nx_g = nx+2*C.num_ghost;
  int ny = C.ny_c;
  int ny_g = ny+2*C.num_ghost;
  // Copy the interior values
  for(int j=0; j < ny; j++)
  {
    for(int i=0; i < nx; i++)
    {
      int xcgint = i+(C.num_ghost+j)*nx_g+C.num_ghost;
      cout << xcgint << endl;
      xc_g[xcgint]= xc[i+j];
      yc_g[xcgint]= yc[i+j];
    }

  }
  for (int i=0; i < nx*ny; i++)
  {
      cout << xc_g[i] << " " << xc_g[i] << endl;
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

  for (int i=0; i < C.nx_c*C.ny_c; i++)
  {
    xc[i] = 0.25*( xn[i] + xn[i+1] + xn[i+nx] + xn[i+nx+1] );
    yc[i] = 0.25*( yn[i] + yn[i+1] + yn[i+nx] + yn[i+nx+1] );
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

void inputMesh(double* xn, double* yn, double* zn, constants C)
{

  // Create the structure constants to contain info needed on all processors

}

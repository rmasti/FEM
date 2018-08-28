//# include "gmsh_io.hpp"
# include "mhdRT.hpp"

int main (int argc, char * argv[])
{
  // Create the structure constants to contain info needed on all processors
  constants C;
  C.f_limiter = 1;
  C.num_ghost = 3;
  C.cfl = 0.5;
  C.nmax = 300;
  C.wint = 50;
  C.pint = 50;
  C.f_mesh = 1;

  int nx_v, ny_v;
  meshSize(&nx_v, &ny_v, C);

  C.nx_c = nx_v-1;
  C.ny_c = ny_v-1;
  double xn[nx_v*ny_v], yn[nx_v*ny_v]; 
  double xc[C.nx_c*C.ny_c], yc[C.nx_c*C.ny_c]; 
  getCoord(xn, yn, xc, yc, C);

  double xc_g[(C.nx_c+2*C.num_ghost)*(C.ny_c+2*C.num_ghost)];
  double yc_g[(C.nx_c+2*C.num_ghost)*(C.ny_c+2*C.num_ghost)];
  extrapCopyCoords(xc_g, yc_g, xc, yc, C);

  cout << nx_v*ny_v << endl;
  return 0;
}

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "gmsh_io.hpp"
# include "mhdRT.hpp"


//int main ( )
int main (int argc, char * argv[])
{
  int *element_node;
  int element_num;
  int element_order;
  string gmsh_filename = "mesh/serial.msh";
  int m;
  int node_num;
  double *node_x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Read data from a file.\n";
  gmsh_size_read ( gmsh_filename, node_num, m, element_num, 
    element_order );
  cout << "\n";
  cout << "  Node data read from file \"" << gmsh_filename << "\"\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Element order = " << element_order << "\n";
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * ) 
    malloc ( element_order * element_num * sizeof ( int ) );
  gmsh_data_read ( gmsh_filename, m, node_num, node_x, 
    element_order, element_num, element_node );
  r8mat_transpose_print_some ( m, node_num, node_x, 
    1, 1, m, 100, "  Coordinates for first 100 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  Connectivity for first 10 elements:" );
//
//  Clean up.
//
  delete [] element_node;
  delete [] node_x;


  return 0;
}

# include "mhdRT.hpp"

int main(int argc, char * argv[])
{
  MPI_Init(& argc, &argv);

  int rank; int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size); 


  MPI_Comm com2d;
  int ni = 3, nj = 5;

  int dim[2] = {4, 3};
  int reorder;
  int period[2], coord[2];
  int id;

  int u, d, l, r;

  
  period[0] = TRUE;
  period[1] = TRUE;
  reorder = FALSE;


  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &com2d); // create 2d comm

/*
  if(rank==7){
  
    printf("\n 0  1  2\n 3  4  5\n 6  7  8\n 9  10  11\n\n");

    MPI_Cart_coords(com2d, rank, 2, coord);

    printf("My rank is %d: My coordinates are %d %d\n", rank, coord[0], coord[1]);

    MPI_Cart_shift(com2d, 0, 1, &u, &d);
    MPI_Cart_shift(com2d, 1, 1, &l, &r);

    printf("My neighbors are left:%d right:%d up:%d down:%d\n", l, r, u ,d);
  
  }

  */



  MPI_Finalize();
  return 0;

}

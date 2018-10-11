# include "mhdRT.hpp"

int main(int argc, char * argv[])
{
  //printf("Run with 12 processors for the example to work\n");
  MPI_Init(& argc, &argv);

  int rank; int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size); 


  double exp=double(log(size))/double(log(2));

  int ni, nj;

  int cnt=size; 
  int n=0;
  /*
   * YANG WAY
   *
   */
  while(cnt >=2)
  {
    cnt = cnt/2;
    n++;
  }
  ni = pow(2, (n/2));
  nj = pow(2, (n/2) + n%2);

  printf("nj = %d, ni = %d, n = %d \n", nj, ni, n);

  /*
  if (exp == ceil(exp))
  {

    if (int(exp) % 2 == 0)
    {
      nj = sqrt(size);
      ni = sqrt(size);

      printf("its even and right: %d\n", int(exp));
      printf("nj = %d, ni = %d\n", nj, ni);
    }
    else
    {
      printf("its odd and right: %d\n", int(exp));
      nj = sqrt(size/2);
      ni = 2*sqrt(size/2);
      printf("nj = %d, ni = %d\n", nj, ni);
    }
  }
  else
  {

    printf("its not right\n");

  }
  */


  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm com2d;

  int dim[2] = {nj, ni};
  int reorder;
  int period[2], coord[2];
  int id;

  int u, d, l, r;


  period[0] = TRUE;//TRUE;
  period[1] = FALSE;//TRUE;
  reorder = FALSE;


  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &com2d); // create 2d comm

  //if(rank==size/2){

  //printf("\n 0  1  2\n 3  4  5\n 6  7  8\n 9  10  11\n\n");

  MPI_Cart_coords(com2d, rank, 2, coord);

  printf("My rank is %d: My coordinates are %d, %d\n", rank, coord[0], coord[1]);

  MPI_Cart_shift(com2d, 0, 1, &u, &d);
  MPI_Cart_shift(com2d, 1, 1, &l, &r);

  printf("My neighbors are left:%d right:%d up:%d down:%d\n", l, r, u ,d);

  //}



  MPI_Finalize();
  return 0;

}

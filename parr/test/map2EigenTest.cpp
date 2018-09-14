#include "mhdRT.hpp"

int main (int argc, char * argv[])
{

  int nrows = 3;
  int ncols = 5;

  Map2Eigen *U = new Map2Eigen(nrows, ncols, NEQ);

  for(int j = 0; j < nrows; j++)
    for(int i = 0; i < ncols; i++)
      for(int eq = 0; eq < NEQ; eq++)
        //  U->Q[eq](j,i) = double((j*ncols+i)*NEQ+eq);
        U->Q_raw[(j*ncols+i)*NEQ+eq] = double((j*ncols+i)*NEQ+eq);
  for(int eq = 0; eq < NEQ; eq++)
  {
    cout << U->Q[eq] << endl;
    cout << endl;
  }

  for(int j = 0; j < nrows; j++)
    for(int i = 0; i < ncols; i++)
      for(int eq = 0; eq < NEQ; eq++)
      {
      
      cout << U->Q_raw[(j*ncols+i)*NEQ+eq] << endl;
      
      cout <<(j*ncols+i)*NEQ+eq<< endl;
      }

  return 0;
}

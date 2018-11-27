#include "mhdRT.hpp"
#include "occa.hpp"

//using namespace occa;

int main(int argc, const char **argv){

  occa::device device("mode: 'CUDA', device_id : 0");

  int nrows = 3;
  int ncols = 5;

  Map2Eigen *U = new Map2Eigen(nrows, ncols, NEQ);
  double* other = new double[nrows*ncols*NEQ];
  for(int j = 0; j < nrows; j++)
    for(int i = 0; i < ncols; i++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q_raw[(j*ncols+i)*NEQ+eq] = double((j*ncols+i)*NEQ+eq);

  cout << U->Q[rhoid] << endl;
 
  occa::memory o_U = device.malloc(nrows*ncols*NEQ*sizeof(double*), U->Q_raw);
  occa::properties props;
  props["defines/NEQ"] = NEQ;

  occa::kernel occaMap2Eigen = occa::buildKernel("test/occaMap2Eigen.okl", "occaMap2Eigen", props);

  occaMap2Eigen(nrows, ncols, o_U);
  
  device.finish();

  o_U.copyTo(U->Q_raw);

  cout << U->Q[rhoid] << endl;

  return 0;
}

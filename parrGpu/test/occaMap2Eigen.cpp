#include "mhdRT.hpp"
#include "occa.hpp"

int main(int argc, const char **argv){

  //occa::setDevice("mode: 'CUDA', "
   //       "device_id : 0");
  //occa::device device("mode: 'Serial'");
  occa::device device("mode: 'CUDA', device_id : 0");

  int nrows = 3;
  int ncols = 5;


  Map2Eigen *U = new Map2Eigen(nrows, ncols, NEQ);
  for(int j = 0; j < nrows; j++)
    for(int i = 0; i < ncols; i++)
      for(int eq = 0; eq < NEQ; eq++)
        U->Q_raw[(j*ncols+i)*NEQ+eq] = double((j*ncols+i)*NEQ+eq);

  cout << U->Q[rhoid] << endl;
 
  occa::memory o_U = device.malloc(nrows*ncols*NEQ*sizeof(double), U->Q_raw);
  //cout << nrows*ncols*NEQ*sizeof(RowMajorMatrixXd) << endl;
  //occa::memory o_U = device.malloc(sizeof(Eigen::MatrixXd)*nrows*ncols, U->Q[rhoid].data());

  occa::properties props;// = occa::deviceProperties();

  props["defines/NEQ"] = NEQ;

  occa::kernel occaMap2Eigen = occa::buildKernel("test/occaMap2Eigen.okl", "occaMap2Eigen", props);

  occaMap2Eigen(nrows, ncols, o_U);

  cout << "\n ________________ MADE IT HERE ____________________\n" << endl;
  
  o_U.copyTo(U->Q_raw);

  device.finish();

  cout << U->Q[rhoid] << endl;

  return 0;
}

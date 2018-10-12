#include "mhdRT.hpp"

int main(int argc, const char **argv){

  //occa::setDevice((std::string) args["options/device"]);
  occa::setDevice("mode: 'Serial'");

  int entries = 5;

  float *a = new float[entries];
  float *b = new float[entries];
  float *ab = new float[entries];

  for(int i = 0; i < entries; i++){
    a[i] = i;
    b[i] = 1-i;
    ab[i] = 0;
  }


  occa::kernel occaAddVector = occa::buildKernel("test/occaAddVector.okl", "occaAddVector");

  occaAddVector(entries, a, b, ab);

  occa::dontSync(b);

  occa::finish();

  for(int i = 0; i < 5; i++)
    cout << i << ": " << ab[i] << endl;

  delete[] a;
  delete[] b;
  delete[] ab;
  return 0;

}

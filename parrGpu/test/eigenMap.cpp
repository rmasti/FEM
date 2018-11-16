#include "mhdRT.hpp"
#include "occa.hpp"

int main(int argc, char * argv[])
{

    int NEQ = 8;
    int ni = 3, nj = 5;
    Map2Eigen *Var = new Map2Eigen(ni, nj, NEQ);
    occa::setDevice("mode: 'CUDA', "
            "device_id : 0");


    for (int i = 0; i < ni*nj*NEQ; i++)
        Var->Q_raw[i] = i;

    Var->Q[rhoid](0,0) = 100;

    occa::setDevice("mode: 'CUDA', "
    occa::setDevice("mode: 'CUDA', "
            "device_id : 0");
            "device_id : 0");
    double *mappoint = (double*) occa::umalloc(ni*nj*NEQ * sizeof(MatrixXd));
    mappoint = Var->Q_raw;

    occa::kernel eigenMap = occa::buildKernel("test/eigenMap.okl", "fillEigenMap");

    eigenMap(Var);

    occa::finish();

    cout << Var->Q[rhoid] << endl;

    return 0;
}

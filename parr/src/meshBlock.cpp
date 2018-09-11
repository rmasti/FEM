#include "mhdRT.hpp"

void meshBlock(const string mesh, const string outputFolder, constants C)
{
  int rank, size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // have the first processor read in data and give it to all other ranks
  if(rank == 0)
  {
    RowMajorMatrixXd xn;
    RowMajorMatrixXd yn;
    RowMajorMatrixXd xc, yc;
    // Read in mesh 
    inputMesh(xn,yn,xc,yc, mesh);

    // extrap to ghost cells
    RowMajorMatrixXd xc_g(xc.rows()+2*C.num_ghost, xc.cols()+2*C.num_ghost);
    RowMajorMatrixXd yc_g(yc.rows()+2*C.num_ghost, yc.cols()+2*C.num_ghost);
    extrapCopyCoords(xc_g, yc_g, xc, yc, C); // get ghost cell coord


    //i or x dir areas the columns of the matrix
    RowMajorMatrixXd Ai(xc.rows(), xn.cols()); // sets (row,col) (j,i)
    //j or y dir areas the rows or the matrix
    RowMajorMatrixXd Aj(xn.rows(), xc.cols()); // sets (row,col) (j,i)
    RowMajorMatrixXd Volume(xc.rows(), xc.cols());

    // SETUP NORMAL VECS
    //Ai normal vector has x and y components (so same size as Ai)
    RowMajorMatrixXd nix(Ai.rows(), Ai.cols());
    RowMajorMatrixXd niy(Ai.rows(), Ai.cols());
    //Aj normal vector has x and y componenets
    RowMajorMatrixXd njx(Aj.rows(), Aj.cols());
    RowMajorMatrixXd njy(Aj.rows(), Aj.cols());

    computeArea(Ai, Aj, xn, yn); // take nodal coords extract interface A // Checked
    computeNormalVectors(nix, niy, njx, njy,
        xn, yn, Ai, Aj); // grab the norm vec 4 computational domain // Checked
    computeVolume(Volume, xn, yn); // Checked

    outputArray(outputFolder, "xc_g", xc_g, 0);
    outputArray(outputFolder, "yc_g", yc_g, 0);
    outputArray(outputFolder, "xc", xc, 0);
    outputArray(outputFolder, "yc", yc, 0);

    outputArray(outputFolder, "Ai", Ai, 0);
    outputArray(outputFolder, "Aj", Aj, 0);
    outputArray(outputFolder, "njx", njx, 0);
    outputArray(outputFolder, "njy", njy, 0);
    outputArray(outputFolder, "nix", nix, 0);
    outputArray(outputFolder, "niy", niy, 0);
    outputArray(outputFolder, "volume", Volume, 0);

    // Now we can send the information that the different blocks need to 
    //cout << xn << endl;


    cout << " Prep to send " << endl;

    int nx = xc.cols();
    int ny = xc.rows();
 
    cout << xn << endl;
    double *xn1d = xn.data();
    //cout << xn.data() << endl;
    //Eigen::Matrix<double, xn.rows(), xn.cols(), Eigen::RowMajor> vertices_;
    //Map<MatrixXd>(xn1d, xn.rows(), xn.cols()) = xn;
    for(int i = 0; i < xn.rows()*xn.cols(); i++)
      cout << xn1d[i] << endl;

    //double ratio = double(nx)/double(ny);
    //cout << nx*ny << " " << nx*ny/size << " "<< sqrt((nx*ny/size)/ratio)  <<" " <<sqrt((nx*ny/size)/ratio)*ratio << endl;
  }
}

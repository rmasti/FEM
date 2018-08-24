#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "mhdRT.hpp"
#define reltol 1.0e-12
using namespace std;


/*
TEST_CASE( "Mtoprim", "[MtoV]" ) {
  constants C;
  C.p0 = 500000.0;
  C.T0 = 2123.9;
  C.gamma = 1.4;
  double M = 1.2343;
  double psi = 1+0.5*(C.gamma-1)*M*M;
  double V1;double V2; double V3; double V4;
  Mtoprim(V1, V2, V3, V4, M, C);
  double pred4 = C.p0*pow(psi,-C.gamma/(C.gamma-1)); 
  double T= C.T0/psi;
  double pred1 = pred4/(R*T);
  double a = compute_soundspeed(C.gamma, pred4, pred1);
  double pred2 = a*M;
  double pred3 = 0.0;
  double tol = mymax(abs(pred1 - V1)/V1, abs(pred2 - V2)/V2);
  tol = mymax(tol,abs(pred4 - V4)/V4) ;
  REQUIRE( tol <= reltol );
}


TEST_CASE( "Sound Speed Calculation", "[a]" ) {
  double p = 658749.9;
  double gamma = 1.4;
  double rho = 87.2990;
  double a = compute_soundspeed(gamma, p, rho);
  double predicted = sqrt(gamma*p/rho);
  REQUIRE( abs(a - predicted)/predicted <= reltol );
}
TEST_CASE( "Test Area function at x = 0.4215", "[A]" ) {
  double A = A_x(0.4215);
  double predicted = 0.50235087913268;
  REQUIRE( abs(A - predicted)/predicted <= reltol );
}

TEST_CASE( "Test M function at x = -0.21315", "[M]" ) {
  double M = M_xinitial(-0.21315);
  double predicted = 0.808165;
  REQUIRE( abs(M - predicted)/predicted <= reltol );
}
*/

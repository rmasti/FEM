/*
 * Main Function file
 * TS: FEM and GPU
 * Written By: Robert Masti
 * 08/28/2018
 */

# include "mhdRT.hpp"


void compute2dFlux(Map2Eigen* F, Map2Eigen* G, const Map2Eigen* U_L , const Map2Eigen* U_R, const Map2Eigen* U_B, const Map2Eigen* U_T, MatrixXd& n_j_xhat, MatrixXd& n_j_yhat, MatrixXd& n_i_xhat, MatrixXd& n_i_yhat, constants C)
{
  int nj = F->Q[rhoid].rows();
  int ni = G->Q[rhoid].cols();
  double ul[NEQ], ur[NEQ], ub[NEQ], ut[NEQ];
  double FFLUX[NEQ];
  double GFLUX[NEQ];

  for(int j = 0; j < nj; j++)
  {
    for(int i = 0; i < ni; i++)
    {
      for(int eq = 0; eq < NEQ; eq++)
      {
        ul[eq] = U_L->Q[eq](j,i);
        ur[eq] = U_R->Q[eq](j,i);
        ub[eq] = U_B->Q[eq](j,i);
        ut[eq] = U_T->Q[eq](j,i);
      }
      computeFlux(FFLUX, ul, ur, n_j_xhat(j,i), n_j_yhat(j,i), 0);
      computeFlux(GFLUX, ub, ut, n_i_xhat(j,i), n_i_yhat(j,i), 1);
      for(int eq = 0; eq < NEQ; eq++)
      {
        F->Q_raw[(j*ni+i)*NEQ+eq] = FFLUX[eq];
        G->Q_raw[(j*ni+i)*NEQ+eq] = GFLUX[eq];
      }
    }
  }
 
  // fix right wall
  for(int j = 0; j < nj; j++)
  {
    for(int eq = 0; eq < NEQ; eq++)
      {
        ul[eq] = U_L->Q[eq](j,ni);
        ur[eq] = U_R->Q[eq](j,ni);

      }
      computeFlux(FFLUX, ul, ur, n_j_xhat(j,ni), n_j_yhat(j,ni), 0);
      for(int eq = 0; eq < NEQ; eq++)
        F->Q[eq](j,ni) = FFLUX[eq];
  }
  // fix left wall

  for(int i = 0; i < ni; i++)
  {
    for(int eq = 0; eq < NEQ; eq++)
      {
        ub[eq] = U_B->Q[eq](nj,i);
        ut[eq] = U_T->Q[eq](nj,i);
 
      }
      computeFlux(GFLUX, ub, ut, n_i_xhat(nj,i), n_i_yhat(nj,i), 1);
      for(int eq = 0; eq < NEQ; eq++)
        G->Q[eq](nj,i) = GFLUX[eq];
  }
}

void computeFlux(double F[], double UA[], double UB[], double& nxhat, double& nyhat, int ForG)
{


  double CA = computeMaxSpeed(UA);
  double CB = computeMaxSpeed(UB);

  double u_A, u_B;

  u_A = (UA[uid]*nxhat + UA[vid]*nyhat)/UA[rhoid]; 
  u_B = (UB[uid]*nxhat + UB[vid]*nyhat)/UB[rhoid]; 

  double lambdaA, lambdaB;
  lambdaA = mymin(u_A, u_B) - mymax(CA, CB);

  lambdaB = mymax(u_A, u_B) + mymax(CA, CB);

  double FA[NEQ], FB[NEQ];

  if (ForG == 0)
  {
    fFlux(FA, UA);
    fFlux(FB, UB);
  }
  else if (ForG == 1)
  {
    gFlux(FA, UA);
    gFlux(FB, UB);
  }
  else
  {
    cerr << "ERROR: Wrong index for direction" << endl;
    exit(-1);
  }

  //for(int eq=0; eq < NEQ; eq++)
   // cout << "Fa = " << FA[eq] << " Fb = " << FB[eq] << endl;
  if (lambdaA > 0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = FA[eq];
  else if (lambdaB < 0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = FB[eq];
  else if (lambdaB >=0 && lambdaA <=0)
    for(int eq = 0; eq < NEQ; eq++)
      F[eq] = (lambdaB*FA[eq] - lambdaA*FB[eq]+lambdaA*lambdaB*(UB[eq]-UA[eq]))/(lambdaB-lambdaA);
  else
  {
    cerr << "ERROR: HLL Flux Calculation Lambda Unphysical" << endl;
    exit(-1);
  }
}
void fFlux(double F[], double U[])
{
  double e;
  double V[NEQ];
  cons2Prim(V, U, NEQ);

  F[rhoid] = V[rhoid]*V[uid]; 
  F[uid] = V[rhoid]*V[uid]*V[uid] - (V[bxid]*V[bxid])/MU0 + V[pid] + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])/(2*MU0);

  F[vid] = V[rhoid]*V[uid]*V[vid]- V[bxid]*V[byid]/MU0;

  F[wid] = V[rhoid]*V[uid]*V[wid] - V[bxid]*V[bzid]/MU0;

  e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

  F[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*V[uid] - (1/MU0)*V[bxid]*(V[uid]*V[bxid]+V[vid]*V[byid]+V[wid]*V[bzid]);

  F[bxid] = 0;

  F[byid] = V[uid]*V[byid]-V[vid]*V[bxid];

  F[bzid] = V[uid]*V[bzid]-V[wid]*V[bxid];
}
void gFlux(double G[], double U[])
{
  double e;
  double V[NEQ];
  cons2Prim(V, U, NEQ);

  G[rhoid] = V[rhoid]*V[vid];

  G[uid] = V[rhoid]*V[uid]*V[vid]- V[bxid]*V[byid]/MU0;

  G[vid] = V[rhoid]*V[vid]*V[vid] - (V[byid]*V[byid])/MU0 + V[pid] + (V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid])/(2*MU0);

  G[wid] = V[rhoid]*V[vid]*V[wid]-V[byid]*V[bzid]/MU0;

  e = V[pid]/(GAMMA-1)+0.5*V[rhoid]*(V[uid]*V[uid]+V[vid]*V[vid]+V[wid]*V[wid]) + 0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]);

  G[pid] = (e+V[pid]+0.5*(1/MU0)*(V[bxid]*V[bxid]+V[byid]*V[byid]+V[bzid]*V[bzid]))*V[vid] - (1/MU0)*V[byid]*(V[uid]*V[bxid]+V[vid]*V[byid]+ V[wid]*V[bzid]);

  G[bxid] = V[vid]*V[bxid]-V[uid]*V[byid];

  G[byid] = 0;

  G[bzid] = V[vid]*V[bzid]-V[wid]*V[byid];
}




void MUSCL(Map2Eigen* U_L, Map2Eigen* U_R, Map2Eigen* U_B, Map2Eigen* U_T,const Map2Eigen* U, constants C)
{
  int ni = U_B->Q[rhoid].cols(); //ncx
  int nj = U_L->Q[rhoid].rows();//ncy

  int ig, jg;
  double theta_L, theta_B, theta_R, theta_T;
  double denom_L, denom_B, denom_R, denom_T;
  double r_L, r_B, r_R, r_T;
  double delta = 1.0e-6;

/*
  printf("\n nj = %d, ni = %d \n", nj, ni);
  printf("\n nj = %ld, ni = %ld \n", U_B->Q[rhoid].rows() , U_B->Q[rhoid].cols());
  printf("\n nj = %ld, ni = %ld \n", U_R->Q[rhoid].rows(), U_R->Q[rhoid].cols());
  printf("\n nj = %ld, ni = %ld \n", U_L->Q[rhoid].rows(), U_L->Q[rhoid].cols());
  printf("\n nj = %ld, ni = %ld \n", U_T->Q[rhoid].rows(), U_T->Q[rhoid].cols());
  printf("\n nj = %ld, ni = %ld \n", U->Q[rhoid].rows(), U->Q[rhoid].cols());
  */
  cout << sizeof U_R->Q_raw << sizeof(double) << endl;
    
  for(int j = 0; j < nj; j++)
  {
  
    for(int i = 0; i < ni; i++)
    {
      ig = i+C.num_ghost;
      jg = j+C.num_ghost;

      for(int eq=0; eq < NEQ; eq++)
      {
        // left right
        
        denom_L = U->Q[eq](jg, ig) - U->Q[eq](jg, ig-1);
        denom_R = U->Q[eq](jg, ig+1) - U->Q[eq](jg, ig-1);
        denom_L = SIGN(1,denom_L)*mymax(delta, abs(denom_L));    
        denom_R = SIGN(1,denom_R)*mymax(delta, abs(denom_R));    
        r_L = (U->Q[eq](jg,ig-1) - U->Q[eq](jg,ig-2))/denom_L;
        r_R = (U->Q[eq](jg,ig) - U->Q[eq](jg,ig-1))/denom_R;
        theta_L = limiter(r_L, C);
        theta_R = limiter(r_R, C);

        denom_B = U->Q[eq](jg, ig) - U->Q[eq](jg-1, ig);
        denom_T = U->Q[eq](jg+1, ig) - U->Q[eq](jg-1, ig);
        denom_B = SIGN(1,denom_B)*mymax(delta, abs(denom_B));    
        denom_T = SIGN(1,denom_T)*mymax(delta, abs(denom_T));    

        r_B = (U->Q[eq](jg-1,ig) - U->Q[eq](jg-2,ig))/denom_B;
        r_T = (U->Q[eq](jg,ig) - U->Q[eq](jg-1,ig))/denom_T;
        theta_B = limiter(r_B, C);
        theta_T = limiter(r_T, C);

        //(U_R->Q_raw[(nj*i+j)*NEQ+eq]) = (U->Q[eq](jg,ig))-0.5*theta_R*((U->Q[eq](jg,ig+1))-(U->Q[eq](jg,ig)));             
        (U_R->Q_raw[(j+i)*NEQ+eq]) = 100.0;

        (U_T->Q_raw[(nj*i+j)*NEQ+eq]) =  (U->Q[eq](jg,ig))+0.5*theta_T*((U->Q[eq](jg+1,ig))-(U->Q[eq](jg,ig)));
        
        (U_L->Q_raw[(nj*i+j)*NEQ+eq]) =  (U->Q[eq](jg,ig-1))+0.5*theta_L*((U->Q[eq](jg,ig))-(U->Q[eq](jg,ig-1)));

        (U_B->Q_raw[(nj*i+j)*NEQ+eq]) =  (U->Q[eq](jg-1,ig))+0.5*theta_B*((U->Q[eq](jg,ig))-(U->Q[eq](jg-1,ig)));

      }
    }
  }
  cout << U_R->Q[rhoid] << endl;

  // add the last layer left right
  for(int j = 0; j < nj; j++)
  {
      ig = ni+C.num_ghost;
      jg = j+C.num_ghost;
    for(int eq = 0; eq < NEQ; eq++)
    {
        denom_L = U->Q[eq](jg, ig) - U->Q[eq](jg, ig-1);
        denom_R = U->Q[eq](jg, ig+1) - U->Q[eq](jg, ig-1);
        denom_L = SIGN(1,denom_L)*mymax(delta, abs(denom_L));    
        denom_R = SIGN(1,denom_R)*mymax(delta, abs(denom_R));    
        r_L = (U->Q[eq](jg,ig-1) - U->Q[eq](jg,ig-2))/denom_L;
        r_R = (U->Q[eq](jg,ig) - U->Q[eq](jg,ig-1))/denom_R;
        theta_L = limiter(r_L, C);
        theta_R = limiter(r_R, C);
    
        U_L->Q[eq](j, ni) =  U->Q[eq](jg,ig-1)+0.5*theta_L*(U->Q[eq](jg,ig)-U->Q[eq](jg,ig-1));
        U_R->Q[eq](j, ni) = U->Q[eq](jg,ig)-0.5*theta_R*(U->Q[eq](jg,ig+1)-U->Q[eq](jg,ig));             
    }
  }
  // add the last layer top bottom
  for(int i = 0; i < ni; i++)
  {
      ig = i+C.num_ghost;
      jg = nj+C.num_ghost;
    for(int eq = 0; eq < NEQ; eq++)
    {
        denom_B = U->Q[eq](jg, ig) - U->Q[eq](jg-1, ig);
        denom_T = U->Q[eq](jg+1, ig) - U->Q[eq](jg-1, ig);
        denom_B = SIGN(1,denom_B)*mymax(delta, abs(denom_B));    
        denom_T = SIGN(1,denom_T)*mymax(delta, abs(denom_T));    
        r_B = (U->Q[eq](jg-1,ig) - U->Q[eq](jg-2,ig))/denom_B;
        r_T = (U->Q[eq](jg,ig) - U->Q[eq](jg-1,ig))/denom_T;
        theta_B = limiter(r_B, C);
        theta_T = limiter(r_T, C);

        U_B->Q[eq](nj, i) =  U->Q[eq](jg-1,ig)+0.5*theta_B*(U->Q[eq](jg,ig)-U->Q[eq](jg-1,ig));
        U_T->Q[eq](nj,i) = U->Q[eq](jg,ig)-0.5*theta_T*(U->Q[eq](jg+1,ig)-U->Q[eq](jg,ig));             
    }
  }
}


double limiter(double& r, constants C)
{
  double theta;
  switch (C.f_limiter)
  {
    case 1:
      theta = 1.0;
    case 2:
      theta = (r+abs(r))/(1+abs(r));
    case 3:
      theta = (r+r*r)/(1+r*r);
    case 4:
      theta = 1.5*(r*r + r) / (1 + r + r*r);
    case 5:
      theta = mymax(0,mymin(2.0*r,mymin(0.5*(1+r),2)));
    case 6:
      theta = mymax(0,mymin(1.0,r));
    case 7:
      theta = mymax(0,mymax(mymin(2.0*r,1),mymin(r,2)));
  }
  return theta;
}            


double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}


double computeTimeStep(
    MatrixXd& Volume,      // input - volume of every cell
    MatrixXd& Ai,          // input - area in i dir 
    MatrixXd& Aj,          // input - area in j dir
    MatrixXd& n_i_xhat,    // input - norm i dir x comp
    MatrixXd& n_i_yhat,    // input - norm i dir y comp
    MatrixXd& n_j_xhat,    // input - norm j dir x comp
    MatrixXd& n_j_yhat,    // input - norm j dir y comp
    Map2Eigen* V,           // input - prim var for speeds
    constants C            // input - constants for the CFL number
    )
{

  int ni = Volume.cols();
  int nj = Volume.rows();

  // need to get the avg normal vectors
  double ni_avg_xhat, ni_avg_yhat;
  double nj_avg_xhat, nj_avg_yhat;
  // get the eigen value which is != MaxSpeed
  double lambda_i, lambda_j;
  double Ai_avg, Aj_avg;

  // add iterables for the cells
  int i_c, j_c;
  int bot = C.num_ghost;// first interior cell in j dir
  int left = C.num_ghost;// first interior cell in i dir
  double dt_min=1e10;
  double maxSpeed;
  double Vmax[NEQ];
  double dt;
  for (int j = 0; j < nj; j++)
  {
    for (int i = 0; i < ni; i++)
    {
      i_c = left + i;
      j_c = bot + j;
      // get the avg norm vecs
      ni_avg_xhat = 0.5*(n_i_xhat(j,i)+n_i_xhat(j,i+1));
      ni_avg_yhat = 0.5*(n_i_yhat(j,i)+n_i_yhat(j,i+1));
      nj_avg_xhat = 0.5*(n_j_xhat(j,i)+n_j_xhat(j+1,i));
      nj_avg_yhat = 0.5*(n_j_yhat(j,i)+n_j_yhat(j+1,i));

      // get avg area
      Ai_avg = 0.5*(Ai(j,i) + Ai(j,i+1));
      Aj_avg = 0.5*(Aj(j,i) + Aj(j+1,i));

      // compute max eigen values
    
      for(int eq = 0; eq < NEQ; eq++)
        Vmax[eq] = V->Q[eq](j_c,i_c);
      
      maxSpeed = computeMaxSpeed(Vmax);
      lambda_i = abs(V->Q[uid](j_c,i_c)*ni_avg_xhat 
          + V->Q[vid](j_c,i_c)*ni_avg_yhat) + maxSpeed;
      lambda_j = abs(V->Q[uid](j_c,i_c)*nj_avg_xhat 
          + V->Q[vid](j_c,i_c)*nj_avg_yhat) + maxSpeed;

      // compute time step non uniform grid
      dt = C.cfl*Volume(j,i) / (lambda_i*Ai_avg + lambda_j*Aj_avg);
      if (dt < dt_min)
        dt_min = dt;
    }
  }
  return dt_min;
}

double computeMaxSpeed(double V[])                                         
{
  double rho, u, v, w, p, bx, by, bz;
  rho = V[rhoid];
  u = V[uid];
  v = V[vid];
  w = V[wid];
  p = V[pid];
  bx = V[bxid];
  by = V[byid];
  bz = V[bzid];
  double cs, ca, cadir;
  double bmag, vmag;
  double maxSpeed;
  bmag = sqrt(bx*bx+by*by+bz*bz);
  vmag = sqrt(u*u+v*v+w*w);
  cs = sqrt(p*GAMMA/rho);
  ca = bmag/sqrt(MU0*rho);
  cadir = ca*((bx*u+by*v+bz*w)/(bmag*vmag+1.0e-12));
  maxSpeed = sqrt(0.5*(ca*ca+cs*cs+sqrt((ca*ca+cs*cs)*(ca*ca+cs*cs)-4*cs*cs*cadir*cadir)));// fast magnetosonic speed

  return maxSpeed;
}                                                                                       



void slipwallBC(
    // This function will allow the normal velocity from the boundary to reflect but simultan
    // eously keep the parrallel to the boundary vel unchanged
    Map2Eigen* V,              // output - Prim var  
    const int Begin[],         // input - Beginning index coord
    const int End[],           // input - Ending index coord 
    const MatrixXd& n_i_xhat,  // input - norm vec i dir x comp
    const MatrixXd& n_i_yhat,  // input - norm vec i dir y comp
    const MatrixXd& n_j_xhat,  // input - norm vec j dir x comp
    const MatrixXd& n_j_yhat,  // input - norm vec j dir y comp
    MatrixXd& T,         // input - temperature at cells
    constants C                // input - constants for num ghosts
    )
{
  int ni_g = V->Q[rhoid].cols();
  int nj_g = V->Q[rhoid].rows();

  int ni = ni_g - 2*C.num_ghost;
  int nj = nj_g - 2*C.num_ghost;

  int i_in_Begin = Begin[1] + C.num_ghost; // 1st interior i dir
  int j_in_Begin = Begin[0] + C.num_ghost; // 1st interior j dir
  int i_in_End = End[1] + C.num_ghost; // shifted by numghost
  int j_in_End = End[0] + C.num_ghost;
  int I, J;
  int sign;
  int BC_index;
  double nx, ny, uvel, vvel;
  double Temp;

  if (i_in_Begin == i_in_End) // Then we are along a vertical boundary loop over j
  {
    if (i_in_Begin == C.num_ghost) // detect whether the ghost cells are on the left (i=3)
    {
      sign = -1; // the face points in the - direction if on the left
      BC_index = Begin[1];
    }
    else if (i_in_Begin == ni + C.num_ghost - 1) // now its on the right wall so sign + and 
    {
      sign = 1;
      BC_index = Begin[1] + 1;
    }
    else
    {
      cerr << "ERROR: Not Applying the Correct SlipWall BC's col!!" << endl;
      exit(1);
    }

    // LOOP OVER THE J-Dir for a vertical boundary!!!!!!!!!!!!!
    for (int j = 0; j <= End[0] - Begin[0]; j++)
    {
      for (int i = 0; i < C.num_ghost; i++)
      {
        I = i_in_Begin + sign*1; // gives the first coln inside the ghost cell
        J = j_in_Begin; // this give first j index
        nx = n_i_xhat(Begin[0]+j, BC_index); // get the normal vec comp in the xdir
        ny = n_i_yhat(Begin[0]+j, BC_index); // get normal vec comp in the y dir
        uvel = V->Q[uid](J+j, i_in_Begin - sign*i); // sign says okay push left for g
        vvel = V->Q[vid](J+j, i_in_Begin - sign*i); // or push right for g
        Temp = T(J+j, i_in_Begin - sign*i); // temperature 

        // get the velocity into comp coords
        V->Q[uid](J+j, I+sign*i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        V->Q[vid](J+j, I+sign*i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);

        // extrapolate P
        V->Q[pid](J+j, I+sign*i) =  2.0*V->Q[pid](J+j, I+sign*(i-1))-V->Q[pid](J+j, I+sign*(i-2));
        V->Q[bxid](J+j, I+sign*i) =  2.0*V->Q[bxid](J+j, I+sign*(i-1))-V->Q[bxid](J+j, I+sign*(i-2));
        V->Q[byid](J+j, I+sign*i) =  2.0*V->Q[byid](J+j, I+sign*(i-1))-V->Q[byid](J+j, I+sign*(i-2));
        V->Q[bzid](J+j, I+sign*i) =  2.0*V->Q[bzid](J+j, I+sign*(i-1))-V->Q[bzid](J+j, I+sign*(i-2));
        V->Q[wid](J+j, I+sign*i) =  2.0*V->Q[wid](J+j, I+sign*(i-1))-V->Q[wid](J+j, I+sign*(i-2));
        T(J+j, I+sign*i) = Temp; // copy temperature extrapolate pressure then compute rho
        V->Q[rhoid](J+j, I+sign*i) = V->Q[pid](J+j, I+sign*i) / (R*T(J+j, I+sign*i));
      }
    }
  }
  else if (j_in_Begin == j_in_End) // now have horizontal boundary
  {
    if (j_in_Begin == C.num_ghost) // at lowere
    {
      sign = -1;
      BC_index = Begin[0];
    }
    else if ( j_in_Begin == nj + C.num_ghost-1 ) // top boundary
    {
      sign = 1;
      BC_index = Begin[0] + 1;
    }
    else
    {
      cerr << "ERROR: Not Applying the Correct SlipWall BC's row!!" << endl;
      exit(1);
    }
    for (int j = 0 ; j < C.num_ghost; j++)
    {
      for (int i = 0 ; i <= End[1]-Begin[1]; i++)
      {
        I = i_in_Begin;          // I is the looping dir
        J = j_in_Begin + sign*1; // j is the num ghost dir
        nx = n_j_xhat(BC_index, Begin[1]+i);
        ny = n_j_yhat(BC_index, Begin[1]+i);
        uvel = V->Q[uid](j_in_Begin - sign*j, I+i); // copy val
        vvel = V->Q[vid](j_in_Begin - sign*j, I+i); // copy val
        Temp = T(j_in_Begin - sign*j, I+i); // copy val
        V->Q[uid](J+sign*j, I+i) = -nx*(uvel*nx+vvel*ny) - ny*(-uvel*ny + vvel*nx);
        V->Q[vid](J+sign*j, I+i) = -ny*(uvel*nx+vvel*ny) + nx*(-uvel*ny + vvel*nx);
        V->Q[pid](J+sign*j, I+i) = 2.0*V->Q[pid](J+sign*(j-1), I+i) - V->Q[pid](J+sign*(j-2), I+i);
        V->Q[wid](J+sign*j, I+i) = 2.0*V->Q[wid](J+sign*(j-1), I+i) - V->Q[wid](J+sign*(j-2), I+i);
        V->Q[bxid](J+sign*j, I+i) = 2.0*V->Q[bxid](J+sign*(j-1), I+i) - V->Q[bxid](J+sign*(j-2), I+i);
        V->Q[byid](J+sign*j, I+i) = 2.0*V->Q[byid](J+sign*(j-1), I+i) - V->Q[byid](J+sign*(j-2), I+i);
        V->Q[bzid](J+sign*j, I+i) = 2.0*V->Q[bzid](J+sign*(j-1), I+i) - V->Q[bzid](J+sign*(j-2), I+i);
        T(J+sign*j, I+i) = Temp;
        V->Q[rhoid](J+sign*j, I+i) = V->Q[pid](J+sign*j, I+i) / (R*T(J+sign*j, I+i));
      }
    }
  }
  else
  {
    cerr << "ERROR: Index Incorrect for SlipWall !!" << endl;
    exit(1);
  }
}

void setBC(
    // Apply BC for cases 2 and 3, God Help Me
    Map2Eigen* V,             // output - prim var
    const MatrixXd& n_i_xhat, // input - norm vec i dir x comp
    const MatrixXd& n_i_yhat, // input - norm vec i dir y comp
    const MatrixXd& n_j_xhat, // input - norm vec j dir x comp
    const MatrixXd& n_j_yhat, // input - norm vec j dir y comp
    MatrixXd& T,        // input - grab the temperature 
    constants C               // input - constants C for case
    )
{
  int ni = V->Q[rhoid].cols()-2*C.num_ghost;
  int nj = V->Q[rhoid].rows()-2*C.num_ghost;

  int Lower_Begin[2]  = {0 , 0};
  int Lower_End[2]    = {0,ni-1};

  int Upper_Begin[2]  = {nj-1, 0};
  int Upper_End[2]    = { nj-1,ni-1};

  int Left_Begin[2]  = {0 , 0};
  int Left_End[2]    = {nj-1,0};

  int Right_Begin[2]  = {0, ni-1};
  int Right_End[2]    = { nj-1,ni-1};

  slipwallBC(V, Upper_Begin, Upper_End, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);
  slipwallBC(V, Lower_Begin, Lower_End, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);

  slipwallBC(V, Left_Begin, Left_End, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);
  slipwallBC(V, Right_Begin, Right_End, n_i_xhat, n_i_yhat, n_j_xhat, n_j_yhat, T, C);

}



void initialize(
    // This function initializes all of the primitive variable and it c
    // hecks for the correct case to initialize with. 
    Map2Eigen* V,           // output - Primitive variables 
    const MatrixXd& xc_g,  // input - xcoord with ghosts
    const MatrixXd& yc_g,  // input - ycoord with ghosts
    constants C            // input - constants for case number
    )
{

  int ni = xc_g.cols(); // n cells in i dir
  int nj = xc_g.rows();  // n cells in j dir
  // check for the same sizes
  if (V->Q[rhoid].cols() != ni || V->Q[rhoid].rows() != nj)
  {
    cerr << "ERROR: Dimension Mismatch in Initialization?!?!" << endl;
    exit(1);
  }
  // create important constants specified from papers
  double x, y, r, thet;
  double rhoL = 1.0;
  double rhoH = 2.0;
  double p0 = 2.5;
  double g = 0.1;

  double vPert = 1.0;
  double rmidd= 0.375;
  double circ = 2.0*PI*rmidd;
  double Bx = 0.0125;
  double lambda = circ/20;
  double randPert;

  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      x = xc_g(j, i);
      y = yc_g(j, i);

      r = sqrt(x*x+y*y);
      thet = atan2(y,x);

      randPert = ((double) rand() / (RAND_MAX));
      if(r < rmidd*(1.0+0.01*cos(10*abs(randPert)*thet/PI + randPert*PI)))
        V->Q[rhoid](j,i) = rhoL;
      else
        V->Q[rhoid](j,i) = rhoH;

      vPert =vPert;
      double decay =exp(-20.0*(r-rmidd)*(r-rmidd)/(0.25*rmidd*rmidd));
      double thetVar = randPert*(2.0*PI/lambda)*thet+randPert*PI;

      V->Q[uid](j,i) = vPert*decay;
      V->Q[vid](j,i) = vPert*decay;
      V->Q[wid](j,i) = 0.0; 
      V->Q[pid](j,i) = p0- V->Q[rhoid](j,i)*g*(r-rmidd);
      V->Q[bxid](j,i) = Bx*cos(thet); 
      V->Q[byid](j,i) = Bx*sin(thet); 
      V->Q[bzid](j,i) = 0.0; 
    }
  }
}




void computeVolume(
    // This function computes the volumes of the interior cell points. 
    MatrixXd& Volume,      // output - Volume of the cells
    const MatrixXd& xn,    // input - x coord for nodal vals
    const MatrixXd& yn     // input - y coord for nodal vals
    )
{
  double dx1, dx2, dy1, dy2; // find the diagonals on the cell
  double cross; // cross product
  for (int i = 0; i < Volume.cols(); i++)
  {
    for (int j = 0; j < Volume.rows(); j++)
    {
      dx1 = xn(j,i) - xn(j+1,i+1); // bot L - top R
      dy1 = yn(j,i) - yn(j+1,i+1); // bot L - top R
      dx2 = xn(j+1,i) - xn(j,i+1); // top L - bot R
      dy2 = yn(j+1,i) - yn(j,i+1); // top L - bot R
      cross = abs(dx1*dy2 - dx2*dy1); // cross product
      Volume(j,i) = 0.5*cross; // assume w=1 volume of cell
    }
  } 
}







void computeNormalVectors(
    // This function takes in the previously calculated areas and
    // Finds the fraction in the x and y physical directions
    MatrixXd& n_i_xhat,    // output - i with A of normal vector in x phys
    MatrixXd& n_i_yhat,    // output - i comp but the yhat phys
    MatrixXd& n_j_xhat,    // output - j with A of normal vector in x phys
    MatrixXd& n_j_yhat,    // output - j dir but with yhat phys
    const MatrixXd& xn,    // input - nodal x coordinate
    const MatrixXd& yn,    // input - nodal y coordinate
    const MatrixXd& Ai,    // input - Area pointed in i (x square grid)
    const MatrixXd& Aj     // input - Area pointed in j (y square grid)
    )
{
  // NOTE i is the columns and j are the rows for A(row,col)
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      // See notes Sec 6 slide 75
      n_i_xhat(j,i) = (yn(j+1,i) - yn(j,i)) / Ai(j,i); 
      n_i_yhat(j,i) = -(xn(j+1,i) - xn(j,i)) / Ai(j,i);
    }
  }
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      n_j_xhat(j,i) = -(yn(j,i+1) - yn(j,i)) / Aj(j,i);
      n_j_yhat(j,i) = (xn(j,i+1) - xn(j,i)) / Aj(j,i);
    }
  }
}



void computeArea(
    // This function computes the area based on the node locations 
    // for both the i facing direction and the j facing direction
    MatrixXd& Ai,          // output - i-direction areas (c, n)size    
    MatrixXd& Aj,          // output - j-direction areas (n, c)size
    const MatrixXd& xn,    // input - x coord for nodes 
    const MatrixXd& yn     // input - y coord for nodes
    )
{
  // loop for Ai
  // NOTE i is the columns and j is the rows
  for (int i = 0; i < Ai.cols(); i++)
  {
    for (int j = 0; j < Ai.rows(); j++)
    {
      double dx = xn(j,i) - xn(j+1,i); // compute bottom - top x
      double dy = yn(j,i) - yn(j+1,i); // compute bottom - top y
      Ai(j,i) = sqrt( dx*dx + dy*dy ); // pythag
    }
  }
  // loop for Aj
  for (int i = 0; i < Aj.cols(); i++)
  {
    for (int j = 0; j < Aj.rows(); j++)
    {
      double dx = xn(j,i) - xn(j,i+1); // compute left - right x
      double dy = yn(j,i) - yn(j,i+1); // compute left - right y
      Aj(j,i) = sqrt( dx*dx + dy*dy);  // pytha
    }
  }
  // NOTE area is the dist times a w into the z dir which is just 1 for 2d
}











void extrapCopyCoords(
    // This function will copy the interior coord of the cells, but
    // will also extrapolate those coords to the ghost cells
    MatrixXd& xc_g,        // output - cell x coord with ghosts 
    MatrixXd& yc_g,        // output - cell y coord with ghosts
    const MatrixXd& xc,    // input - cell x coord without ghosts
    const MatrixXd& yc,    // input - cell y coord without ghosts
    constants C            // input - constants C for num ghost
    )
{
  int ni = xc.cols();
  int nj = xc.rows();
  // copy the interior values for both x and y
  for (int i = 0; i < ni; i++)
  {
    for (int j = 0; j < nj; j++)
    {
      int i_g = i+C.num_ghost; // set iter to 3+i
      int j_g = j+C.num_ghost; // set iter to 3+j
      xc_g(j_g, i_g) = xc(j,i);
      yc_g(j_g, i_g) = yc(j,i);
    }
  }
  int bot = C.num_ghost - 1; // 1st ghost cell based off 
  int top = C.num_ghost + nj; // top index start
  int left = C.num_ghost -1; // same as bot
  int right = C.num_ghost + ni; // 1st ghost on right
  // extrapolate to the ghost cells now

  // j dir which means the interior need extrap
  for (int i = C.num_ghost; i < ni+C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bottom index - j is the filling of ghost
      xc_g(bot-j,i) = 2.0*xc_g(bot-j+1,i) - xc_g(bot-j+2,i);
      yc_g(bot-j,i) = 2.0*yc_g(bot-j+1,i) - yc_g(bot-j+2,i);
      // top index top + j is filling ghost
      xc_g(top+j,i) = 2.0*xc_g(top+j-1,i) - xc_g(top+j-2,i);
      yc_g(top+j,i) = 2.0*yc_g(top+j-1,i) - yc_g(top+j-2,i);
    }
  }
  // i dir which means the interior need extrap
  for (int i = 0; i < C.num_ghost; i ++)
  {
    for (int j = C.num_ghost; j < C.num_ghost+nj; j++)
    {
      // left index - i is the filling of ghost
      xc_g(j,left-i) = 2.0*xc_g(j,left-i+1) - xc_g(j,left-i+2);
      yc_g(j,left-i) = 2.0*yc_g(j,left-i+1) - yc_g(j,left-i+2);
      // right index + i is the filling of ghost
      xc_g(j,right+i) = 2.0*xc_g(j,right+i-1) - xc_g(j,right+i-2);
      yc_g(j,right+i) = 2.0*yc_g(j,right+i-1) - yc_g(j,right+i-2);
    }
  }
  //set to 1 out the corners for sake of clarity.
  for (int i = 0; i < C.num_ghost; i++)
  {
    for (int j = 0; j < C.num_ghost; j++)
    {
      // bot left
      xc_g(j,i) = 1.0; yc_g(j,i) = 1.0;
      // bot right
      xc_g(j,right+i) = 1.0; yc_g(j,right+i) = 1.0;
      // top right
      xc_g(top+j,right+i) = 1.0; yc_g(top+j,right+i) = 1.0;
      // top left
      xc_g(top+j,i) = 1.0; yc_g(top+j,i) = 1.0;
    }
  }
}

void inputMesh(
    MatrixXd& xn,         // Output - x coord nodes
    MatrixXd& yn,         // Output - y coord nodes
    MatrixXd& xc,         // Output - x coord centers
    MatrixXd& yc,         // Output - y coord centers
    const string mesh)    // input - mesh location
{
  ifstream infile;
  infile.open(mesh.c_str());

  if (!infile){
    cerr << "ERROR: Unable to open the mesh file" << endl;
    exit(1);
  }

  double readval;
  string line;
  int vals[3];
  for(int i = 0; i < 3; i++ )
  {
    getline(infile, line);
    istringstream streamA(line);
    streamA >> readval;
    vals[i] = readval;
  }
  int ntot = vals[0];
  int nx = vals[1];
  int ny = vals[2];

  double dat[ntot*2];
  int i = 0;
  while (infile.good())
    while (getline(infile, line))
    {
      istringstream streamA(line);
      while(streamA >> readval)
        dat[i] = readval;
      i++;
    } 
  xn.resize(ny, nx); //number of rows, number of columns
  yn.resize(ny, nx); //number of rows, number of columns
  xc.resize(ny-1, nx-1);
  yc.resize(ny-1, nx-1);

  int offset = ntot*2;
  for(int j = 0; j < ny; j++)
    for(int i = 0; i < nx; i++)
    {
      int xInd = j*nx+i;
      int yInd = j*nx+i+nx*ny;
      xn(j, i) = dat[xInd];
      yn(j, i) = dat[yInd];
    }

  for (int j = 0; j < xc.rows(); j++)
    for (int i = 0 ; i < xc.cols(); i++)
    {
      // take average of all four xvals from corn
      // order: botL + topL + botR + topR 
      xc(j,i) = 0.25*(xn(j,i) + xn(j+1,i) + 
          xn(j,i+1) + xn(j+1,i+1));
      yc(j,i) = 0.25*(yn(j,i) + yn(j+1,i) + 
          yn(j,i+1) + yn(j+1,i+1));
    }
}

void outputArrayMap(
    // This function will output any eigen matrix into a file
    string Address,        // input - folder location
    string FileName,       // input - File
    const Map2Eigen* out,         //input - matrix
    int n                  //input - iteration
    )
{
  string outString = FileName + "_rho";
  MatrixXd outSingle = out->Q[rhoid];

  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[uid]; outString = FileName+"_u" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[vid]; outString = FileName+"_v" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[pid]; outString = FileName+"_p" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[bxid]; outString = FileName+"_bx" ;
  outputArray(Address, outString, outSingle, n);
  outSingle = out->Q[byid]; outString = FileName+"_by" ;
  outputArray(Address, outString, outSingle, n);
}

void outputArray(
    // This function will output any eigen matrix into a file
    string Address,        // input - folder location
    string FileName,       // input - File
    MatrixXd& out,         //input - matrix
    int n                  //input - iteration
    )
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Suffix = ".txt";
  string Sub = "-";
  Address = Address + "/" + FileName + Sub + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component

  // Alot of options for outputting this is just nice to see 
  // where lines end and start
  //IOFormat CleanFmt(14,0,", ", "\n", "[", "]");
  IOFormat CleanFmt(14);
  outfile << out.transpose().format(CleanFmt) << endl; // output trans
  //outfile << setprecision(14) << out << endl; // output trans
}

void primToCons(
    Map2Eigen* U,           // output - Conserved vars
    const Map2Eigen* V)         // input - prim vars
{
  int ni = U->Q[rhoid].cols();
  int nj = U->Q[rhoid].rows();
  for(int j = 0; j < nj; j++)
    for(int i = 0; i < ni; i++)
    {
      U->Q_raw[(j*ni+i)*NEQ+rhoid] = V->Q[rhoid](j,i); //U1 = V1
      U->Q_raw[(j*ni+i)*NEQ+uid] = V->Q[rhoid](j,i)*(V->Q[uid](j,i));
      U->Q_raw[(j*ni+i)*NEQ+vid] = V->Q[rhoid](j,i)*(V->Q[vid](j,i));
      U->Q_raw[(j*ni+i)*NEQ+wid]= V->Q[rhoid](j,i)*(V->Q[wid](j,i));
      U->Q_raw[(j*ni+i)*NEQ+pid] = V->Q[pid](j,i)/(GAMMA - 1.0) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[uid](j,i))*(V->Q[uid](j,i))) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[vid](j,i))*(V->Q[vid](j,i))) 
        + 0.5*V->Q[rhoid](j,i)*((V->Q[wid](j,i))*(V->Q[wid](j,i))) 
        + 0.5*V->Q[bxid](j,i)*((V->Q[bxid](j,i)))/MU0 
        + 0.5*V->Q[byid](j,i)*(V->Q[byid](j,i))/MU0 
        + 0.5*V->Q[bzid](j,i)*(V->Q[bzid](j,i))/MU0;
    
      U->Q_raw[(j*ni+i)*NEQ+bxid] = V->Q[bxid](j,i); 
      U->Q_raw[(j*ni+i)*NEQ+byid] = V->Q[byid](j,i); 
      U->Q_raw[(j*ni+i)*NEQ+bzid] = V->Q[bzid](j,i); 
    }
}



void consToPrim(
    Map2Eigen* V,           // output - Conserved vars
    const Map2Eigen* U)        // input - prim vars
{
  int ni = V->Q[rhoid].cols();
  int nj = V->Q[rhoid].rows();
  for(int j = 0; j < nj; j++)
    for(int i = 0; i < ni; i++)
    {
      V->Q_raw[(j*ni+i)*NEQ+rhoid] = U->Q[rhoid](j,i); //U1 = V1
      V->Q_raw[(j*ni+i)*NEQ+uid] = U->Q[uid](j,i)/(U->Q[rhoid](j,i));
      V->Q_raw[(j*ni+i)*NEQ+vid] = U->Q[vid](j,i)/(U->Q[rhoid](j,i));
      V->Q_raw[(j*ni+i)*NEQ+wid] = U->Q[wid](j,i)/(U->Q[rhoid](j,i));
      V->Q_raw[(j*ni+i)*NEQ+pid] = (GAMMA - 1.0) *(U->Q[pid](j,i)  
        - 0.5*U->Q[uid](j,i)*((U->Q[uid](j,i))/(U->Q[rhoid](j,i))) 
        - 0.5*U->Q[vid](j,i)*((U->Q[vid](j,i))/(U->Q[rhoid](j,i))) 
        - 0.5*U->Q[wid](j,i)*((U->Q[wid](j,i))/(U->Q[rhoid](j,i))) 
        - 0.5*U->Q[bxid](j,i)*(U->Q[bxid](j,i))/MU0 
        - 0.5*U->Q[byid](j,i)*(U->Q[byid](j,i))/MU0 
        - 0.5*U->Q[bzid](j,i)*(U->Q[bzid](j,i))/MU0);
    
      V->Q_raw[(j*ni+i)*NEQ+bxid] = U->Q[bxid](j,i); 
      V->Q_raw[(j*ni+i)*NEQ+byid] = U->Q[byid](j,i); 
      V->Q_raw[(j*ni+i)*NEQ+bzid] = U->Q[bzid](j,i); 
    }
}

void cons2Prim(
    double V[],           // Output: Primitive variables
    double U[],           // Input: Conservative variables
    int size)             // Input: size of array must be evenly dividable by NEQ
{
  int index; 
  double U0[NEQ];
  for(int i = 0; i < size/NEQ; i++)
  {
    index = i*NEQ;
    for(int eq = 0; eq < NEQ; eq++)
      U0[eq] = U[index+eq], V[index+eq] = U0[eq];

    V[index+uid] = U0[uid]/U0[rhoid]; //rhou/rho
    V[index+vid] = U0[vid]/U0[rhoid]; //rhou/rho
    V[index+wid] = U0[wid]/U0[rhoid]; //rhou/rho
    V[index+pid] = (GAMMA-1.0)*(U0[pid] - 0.5*(U0[uid]*U0[uid] + U0[vid]*U0[vid] + U0[wid]*U0[wid])/U0[rhoid] - 0.5*((U0[bxid]*U0[bxid] + U0[byid]*U0[byid] + U0[bzid]*U0[bzid]))/MU0);
  }
}



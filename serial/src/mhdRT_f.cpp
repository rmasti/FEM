# include "mhdRT.hpp"

void compute2dFlux(double F[], double G[], double Ul[], double Ur[], double Ub[], double Ut[], double njx[], double njy[], double nix[], double niy[], constants C)
{
  double UL_L[NEQ], UR_L[NEQ], UB_B[NEQ], UT_B[NEQ];
  double UL_R[NEQ], UR_R[NEQ], UB_T[NEQ], UT_T[NEQ];
  double FLUX[NEQ];
  int nx = C.nx_c+1;
  int ri, li, bi, ti;
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int c = i*C.nx_c + j;
      li = i*(nx)+j;
      ri = i*(nx)+j+1;
      ti = (i+1)*(C.nx_c)+j;
      bi = (i)*C.nx_c+j;

      for(int eq = 0; eq < NEQ; eq++)
      {
        // left and bottom faces
        UL_L[eq] = Ul[li*NEQ+eq]; UR_L[eq] = Ur[li*NEQ+eq];
        UB_B[eq] = Ub[bi*NEQ+eq]; UT_B[eq] = Ut[bi*NEQ+eq];

        // right and top faces
        UL_R[eq] = Ul[ri*NEQ+eq]; UR_R[eq] = Ur[ri*NEQ+eq];
        UB_T[eq] = Ub[ti*NEQ+eq]; UT_T[eq] = Ut[ti*NEQ+eq];
      }
      cout << endl;
      //printf("c = %d, li = %d, ri = %d, bi = %d, ti = %d\n", c, li, ri, bi, ti);

      computeFlux(FLUX, UL_L, UR_R, njx[li], njy[li], 0);

      //cout << FLUX[0] << " " << endl;


      /*
         CR = computeMaxSpeed(UR); CL = computeMaxSpeed(UL);
         CB = computeMaxSpeed(UB); CL = computeMaxSpeed(UT);
         u_L = (UL[cneq+uid]*njx[li] + UL[cneq+vid]*njy[li])/UL[cneq+rhoid]; 
         u_R = (UR[cneq+uid]*njx[li] + UR[cneq+vid]*njy[li])/UR[cneq+rhoid]; 
         */

    }
  }
}

void computeFlux(double F[], double UA[], double UB[], double& nxhat, double& nyhat, int ForG)
{
  for(int eq = 0; eq < NEQ ; eq++)
  {
    //printf("UA[%d] = %lf, ", eq, UA[eq]);


  }
  //cout << endl;
  double CA = computeMaxSpeed(UA);
  double CB = computeMaxSpeed(UB);
  // cout << CA << " " << CB << endl;

}


void computeSource(double S[], double U[], double thetc[], constants C)
{
  // rt source term rho u rho v rho uu rho vv
  double g = C.g;
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {    
      int c = i*C.nx_c + j;
      int cneq = c*NEQ;
      S[cneq+uid] = U[cneq+rhoid]*g*cos(thetc[c]);//rho g costhet u src 
      S[cneq+vid] = U[cneq+rhoid]*g*sin(thetc[c]);//rho g sinthet v src
      S[cneq+pid]=g*(U[cneq+uid]*cos(thetc[c])+U[cneq+vid]*sin(thetc[c])); 
    }
  }
}

void MUSCL(double lout[], double rout[], double bout[], double tout[], double Ul_g[], double Ur_g[], double Ub_g[], double Ut_g[], double U[], constants C)
{

  int nx = C.nx_c + 1;
  int ny = C.ny_c + 1;
  double ULR[NEQ*5]; //5 point stencil
  double UBT[NEQ*5]; //5 point stencil

  double UO[NEQ]; //5 point stencil

  double theta_L, theta_R, theta_B, theta_T;
  double delta = 1.0e-6;
  int li, ri;
  int bi, ti;

  // left right top bottom
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int g = i*C.num_ghost;
      int gv = j*C.num_ghost;
      int c = i*C.nx_c + j;
      li = i*(nx)+j;
      ri = i*(nx)+j+1;
      ti = (i+1)*(C.nx_c)+j;
      bi = (i)*C.nx_c+j;

      //printf("c = %d, li = %d, ri = %d, bi = %d, ti = %d\n", c, li, ri, bi, ti);
      for(int k = 0; k < 5; k++)
      {
        for(int eq = 0; eq < NEQ; eq++)
        {
          ULR[k*NEQ+eq] = U[(c+k-2)*NEQ+eq]; // get initial lr
          UBT[k*NEQ+eq] = U[(c+(k-2)*C.nx_c)*NEQ+eq]; // get initial tb
        }
      }

      // check if cell is bot or top side
      // THIS WORKS 
      if(i < 2) // on the bottom
        if(i == 0) // on the very bottom
          for(int k = 1; k >= 0; k--)
            for(int eq = 0; eq < NEQ; eq++)
              UBT[k*NEQ+eq] = Ub_g[(gv+(1-k))*NEQ+eq];
        else
          for(int eq = 0; eq < NEQ; eq++)
          {
            int k = 0;
            UBT[k*NEQ+eq] = Ub_g[(gv+k)*NEQ+eq];
          }

      if(j < 2) // on the left
        if(j == 0) // on left side 1st cell
          for(int k = 1; k >= 0; k--)
            for(int eq = 0; eq < NEQ; eq++)
              ULR[k*NEQ+eq] = Ul_g[(g+(1-k))*NEQ+eq];
        else       // on left side second cell
          for(int eq = 0; eq < NEQ; eq++)
          {
            int k = 0; 
            ULR[k*NEQ+eq] = Ul_g[(g+k)*NEQ+eq];
          }

      if( i >= C.ny_c-2) // top side
        if( i == C.ny_c-2) // top side 2 layers in
        {
          int k = 4;
          for(int eq=0; eq < NEQ; eq++)
            UBT[k*NEQ+eq] = Ut_g[eq+NEQ*(gv+k-4)];
        }
        else
          for(int k = 3; k < 5; k++)
            for(int eq = 0; eq < NEQ; eq++)
              UBT[k*NEQ+eq] = Ut_g[eq+NEQ*(gv+k-3)];

      if( j >= C.nx_c-2) // on right side
        if(j-(C.nx_c-1) < 0)// on right side 2 cells in 
        {
          int k = 4;
          for(int eq=0; eq < NEQ; eq++)
            ULR[k*NEQ+eq] = Ur_g[eq+NEQ*(g+k-4)];
        }
        else // on right side 1 cell in
          for(int k = 3; k < 5; k++)
            for(int eq = 0; eq < NEQ; eq++)
              ULR[k*NEQ+eq] = Ur_g[eq+NEQ*(g+k-3)];

         outputArray("output", "ULR", ULR, (sizeof ULR/sizeof *ULR), c);
         /*
         outputArray("output", "UBT", UBT, (sizeof UBT/sizeof *UBT), c);
         */

      //outputArray("output", "ULR", ULR, (sizeof ULR/sizeof *ULR), c);
      for(int eq = 0; eq < NEQ; eq++)
      {

        int fi = 2; // left and bottom int 

        getThetaExtrap(theta_L, theta_R, ULR, fi, eq, C);
        getThetaExtrap(theta_B, theta_T, UBT, fi, eq, C);

        /*
        lout[li*NEQ+eq]= ULR[NEQ*(fi-1)+eq] + 0.5*theta_L*(ULR[NEQ*(fi)+eq] - ULR[NEQ*(fi-1)+eq]);

        rout[li*NEQ+eq] = ULR[NEQ*(fi)+eq] - 0.5*theta_R*(ULR[NEQ*(fi+1)+eq] - ULR[NEQ*(fi)+eq]);

        bout[bi*NEQ+eq] = UBT[NEQ*(fi-1)+eq] + 0.5*theta_B*(UBT[NEQ*(fi)+eq] - UBT[NEQ*(fi-1)+eq]);
        tout[bi*NEQ+eq] = UBT[NEQ*(fi)+eq] - 0.5*theta_T*(UBT[NEQ*(fi+1)+eq] - UBT[NEQ*(fi)+eq]);

        fi = 3; // right and top int

        getThetaExtrap(theta_L, theta_R, ULR, fi, eq, C);
        getThetaExtrap(theta_B, theta_T, UBT, fi, eq, C);

        lout[ri*NEQ+eq] = ULR[NEQ*(fi-1)+eq] + 0.5*theta_L*(ULR[NEQ*(fi)+eq] - ULR[NEQ*(fi-1)+eq]);
        rout[ri*NEQ+eq] = ULR[NEQ*(fi)+eq] - 0.5*theta_R*(ULR[NEQ*(fi+1)+eq] - ULR[NEQ*(fi)+eq]);

        bout[ti*NEQ+eq] = UBT[NEQ*(fi-1)+eq] + 0.5*theta_B*(UBT[NEQ*(fi)+eq] - UBT[NEQ*(fi-1)+eq]);
        tout[ti*NEQ+eq] = UBT[NEQ*(fi)+eq] - 0.5*theta_T*(UBT[NEQ*(fi+1)+eq] - UBT[NEQ*(fi)+eq]);

        */
        //printf("c = %d, li = %d, ri = %d, bi = %d, ti = %d\n", c, li, ri, bi, ti);
      }
    }
  }
}

void getThetaExtrap(double& theta_A, double& theta_B, double U5[], int& ind, int& eq, constants& C)
{
  double denom_L, denom_R;
  double r_L, r_R;
  double delta = 1.0e-6;

  denom_L = (U5[NEQ*(ind)+eq]-U5[NEQ*(ind-1)+eq]);
  denom_R = (U5[NEQ*(ind+1)]-U5[NEQ*(ind)+eq]);
  denom_L = SIGN(1,denom_L)*mymax(delta, abs(denom_L));
  denom_R = SIGN(1,denom_R)*mymax(delta, abs(denom_R));
  r_L = (U5[NEQ*(ind-1)+eq] - U5[NEQ*(ind-2)+eq])/denom_L;
  r_R = (U5[NEQ*(ind)+eq] - U5[NEQ*(ind-1)+eq])/denom_R;
  theta_A = limiter(r_L, C);
  theta_B = limiter(r_R, C);
}

double limiter(double& r, constants C)
{
  double theta;
  /*
     %   theta = 1 no limiter
     %   theta = 2 van leer
     %   theta = 3 van albaada
     %   theta = 4 ospre
     %   theta = 5 monotized central least diffusive
     %   theta = 6 MINMOD
     %   theta = 7 superbee
     */

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

double computeTimeStep(double volume[], double Aj[], double Ai[], double njx[], double njy[], double nix[], double niy[], double V[], constants C)
{
  int nx = C.nx_c + 1;
  int ny = C.ny_c + 1;
  double dto=1e8;
  double dt=1e8;

  double nixhat, niyhat; // avg vals
  double njxhat, njyhat; // avg vals
  // to get eigen val
  double lambda_j, lambda_i;
  double Aj_avg, Ai_avg;

  int c;
  int li, ri;
  int ti, bi;

  double V0[NEQ];
  double maxSpeed;

  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++ )
    {
      c = i*C.nx_c + j;
      li = i*nx + j;
      ri = li + 1;
      bi = i*C.nx_c + j;
      ti = (i+1)*C.nx_c + j;
      // left right avg
      njxhat = 0.5*(njx[li] + njx[ri]);
      njyhat = 0.5*(njy[li] + njy[ri]);
      nixhat = 0.5*(nix[bi] + nix[ti]);
      niyhat = 0.5*(niy[bi] + niy[ti]);

      // A avg
      Aj_avg = 0.5*(Aj[li]+Aj[ri]);
      Ai_avg = 0.5*(Ai[bi]+Ai[ti]);

      for(int eq = 0; eq < NEQ; eq++)
        V0[eq] = V[c*NEQ+eq]; 

      maxSpeed = computeMaxSpeed(V0); 
      lambda_j = abs(V0[uid]*njxhat + V0[vid]*njyhat) + maxSpeed;
      lambda_i = abs(V0[uid]*nixhat + V0[vid]*niyhat) + maxSpeed;


      dt = C.cfl*volume[c]/(lambda_j*Aj_avg + lambda_i*Ai_avg);

      if (dt < dto)
        dto = dt;
      //printf("cell = %d, l = %d, r = %d, b = %d, t = %d \n", c, li, ri, bi, ti);
    }
  }
  dt = dto;

  return dt;
}

void setBC(double Vl_g[], double Vr_g[], double Vt_g[], double Vb_g[], double njx[],double njy[], double nix[], double niy[], double V[], constants C)
{
  int nx = C.nx_c + 1;
  int ny = C.ny_c + 1;

  switch (C.f_case)
  {
    case 1: // periodic left and right
      double V0[NEQ];
      for(int i = 0; i < C.ny_c; i++)
      {
        for(int j = 0; j < C.num_ghost; j++)
        {
          int lInt = i*C.nx_c + j;
          int rInt = i*(C.nx_c) + (C.nx_c-1) - j;
          int g = i*C.num_ghost+j;

          //printf("lInt = %d, rInt = %d, g = %d\n", lInt, rInt, g);

          for(int eq = 0; eq < NEQ; eq++)
          {
            double coeff;
            // Must flip u and v because of the geometry cyl coord half circ only
            if ( eq == uid || eq==vid || eq==bxid || eq==byid)
              coeff = -1.0;
            else
              coeff = 1.0;

            // flips the things needed to be flipped
            Vl_g[g*NEQ+eq] = coeff*V[rInt*NEQ+eq]; // left grabs right
            Vr_g[g*NEQ+eq] = coeff*V[lInt*NEQ+eq]; // right grabs left
          }
        }
      }
      // Slip wall top and bottom see leveque page 486-487
      for(int j=0; j < C.nx_c; j++)
      {
        int tFnt = C.nx_c*(ny-1) + j;
        int bFnt =  j ;
        double txhat = nix[tFnt];
        double bxhat = nix[bFnt];
        double tyhat = niy[tFnt];
        double byhat = niy[bFnt];

        for(int i=0; i < C.num_ghost; i++)
        {
          int g = C.num_ghost*j + i;
          int bInt = (C.nx_c*i) + j ;
          int tInt = C.nx_c*(C.ny_c-i-1) + j;

          for(int eq = 0; eq < NEQ; eq++)
          {
            Vt_g[g*NEQ+eq] = V[tInt*NEQ+eq]; // top grabs top
            Vb_g[g*NEQ+eq] = V[bInt*NEQ+eq]; // bot grabs bot
          }
          // Apply the slipwall with perfect conduction
          double u = V[tInt*NEQ+uid];
          double v = V[tInt*NEQ+vid];
          double bx = V[tInt*NEQ+bxid];
          double by = V[tInt*NEQ+byid];

          // No vel penetration
          Vt_g[g*NEQ+uid] = -txhat*(u*txhat+v*tyhat) - tyhat*(-u*tyhat + v*txhat);
          Vt_g[g*NEQ+vid] = -tyhat*(u*txhat+v*tyhat) + txhat*(-u*tyhat + v*txhat);

          Vb_g[g*NEQ+uid] = -bxhat*(u*bxhat+v*byhat) - byhat*(-u*byhat + v*bxhat);
          Vb_g[g*NEQ+vid] = -byhat*(u*bxhat+v*byhat) + bxhat*(-u*byhat + v*bxhat); 

          // no b penetration
          Vt_g[g*NEQ+bxid] = -txhat*(bx*txhat+by*tyhat) - tyhat*(-bx*tyhat + by*txhat);
          Vt_g[g*NEQ+byid] = -tyhat*(bx*txhat+by*tyhat) + txhat*(-bx*tyhat + by*txhat);

          Vb_g[g*NEQ+bxid] = -bxhat*(bx*bxhat+by*byhat) - byhat*(-bx*byhat + by*bxhat); 
          Vb_g[g*NEQ+byid] = -byhat*(bx*bxhat+by*byhat) + bxhat*(-bx*byhat + by*bxhat); 
        }
      }

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

void prim2Cons(
    double U[],           // Output: Conserved variables
    double V[],           // Input: Primitive variables
    int size)             // Input: size of array must be evenly dividable by NEQ
{
  int index;
  double V0[NEQ];
  for(int i = 0; i < size/NEQ; i++)
  {
    index = i*NEQ;
    for(int eq = 0; eq < NEQ; eq++)
      V0[eq] = V[index+eq], U[index+eq] = V0[eq];

    U[index+uid] = V0[rhoid]*V0[uid]; // rhou
    U[index+vid] = V0[rhoid]*V0[vid]; // rhov
    U[index+wid] = V0[rhoid]*V0[wid]; // rhow
    U[index+pid] = V0[pid]/(GAMMA-1.0) + 0.5*V0[rhoid]*(V0[uid]*V0[uid] + V0[vid]*V0[vid] + V0[wid]*V0[wid]) + 0.5*((V0[bxid]*V0[bxid] + V0[byid]*V0[byid] + V0[bzid]*V0[bzid]))/MU0; // e 
  }
}

void initialize(
    double V[],           // Output: computes the volume of each sell doing cross product
    double rc[],          // Input: radius coordinator
    double thetc[],       // Input: theta coordinator
    constants C)          // Input: constants for cell sizes, initial cases
{
  double V0[NEQ];
  double r, thet;
  int cIndex;
  switch (C.f_case)
  {
    case 1:
      double rhoL = 1.0;
      double rhoH = 2.0;
      double p0 = 2.5;
      double g = C.g;
      double vPert = 0.01;
      double bparr0 = 0.0*sqrt(4.0*PI);
      double circ = 2.0*PI*0.375;
      double lambda = circ/20;
      double randPert;
      V0[wid] = 0.0;

      for(int i = 0; i < C.ny_c; i++)
      {
        for(int j = 0; j < C.nx_c; j++)
        {
          cIndex = i*C.nx_c + j;
          r = rc[cIndex];
          thet = thetc[cIndex];
          if(r < 0.375)
          {
            V0[rhoid] = rhoL;
            V0[pid] = p0-V0[rhoid]*g*(r-0.375);
          }
          else
          {
            V0[rhoid] = rhoH;
            V0[pid] = p0 - V[rhoid]*g*(r-0.375);
          }; 
          randPert = ((double) rand() / (RAND_MAX));
          double decay =exp(-20.0*(r-0.375)*(r-0.375)/(0.125*0.125)); 
          double thetVar = randPert*(2.0*PI/lambda)*thet+randPert*PI;
          V0[uid] = vPert*cos(thetVar)*decay; 
          V0[vid] = vPert*sin(thetVar)*decay; 
          V0[bxid] = bparr0*cos(thet);
          V0[byid] = bparr0*sin(thet);
          V0[bzid] = 0.0;

          // Fill in V
          for(int eq = 0; eq < NEQ; eq++)
            V[cIndex*NEQ+eq] = V0[eq];
        }
      }
  }
}

void computeVolume(
    double volume[],      // Output: computes the volume of each sell doing cross product
    double xn[],          // Input: node coords x dir
    double yn[],          // Input: node coords y dir
    constants C)          // Input: constants for cell sizes
{
  double dx1, dx2, dy1, dy2, cross;
  int nx = C.nx_c + 1;
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int cIndex = i*C.nx_c+j;
      int nIndex = i*nx+j;

      //bot L - top R
      dx1 = xn[nIndex] - xn[nIndex+nx+1];
      dy1 = yn[nIndex] - yn[nIndex+nx+1];
      //top L - bot R
      dx2 = xn[nIndex+nx] - xn[nIndex+1];
      dy2 = yn[nIndex+nx] - yn[nIndex+1];

      cross = abs(dx1*dy2 - dx2*dy1);
      volume[cIndex] = 0.5*cross;
    }
  }
}

void computeAreasAndNormalVectors(
    double njx[],         // Output: normal vector pointing in j increase col dir xhat
    double njy[],         // Output: normal vector pointing in j col dir yhat
    double nix[],         // Output: normal vector pointing in i row dir xhat
    double niy[],         // Output: normal vector pointing in i row dir yhat
    double Aj[],          // Output: Area for j col facing dir nc x ni
    double Ai[],          // Output: Area for i col facing dir ni x nc
    double xn[],          // Input: x node coords
    double yn[],          // Input: y node coords
    constants C)          // Input: Constants for cell sizes
{
  int nx = C.nx_c+1; // number of interfaces/nodes x dir
  int ny = C.ny_c+1; // and y dir

  // loop for Aj or x dir
  for(int i = 0; i < C.ny_c; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      int nIndex = i*nx+j;
      double dx = xn[nIndex] - xn[nIndex+nx];
      double dy = yn[nIndex] - yn[nIndex+nx];
      Aj[nIndex] = sqrt(dx*dx + dy*dy);
      njx[nIndex] = dy / Aj[nIndex];
      njy[nIndex] = -1.0*dx / Aj[nIndex];
    }
  }
  // loop for Ai or the ydir
  for(int i = 0; i < ny; i++)
  {
    for(int j = 0; j < C.nx_c; j++)
    {
      int nIndex = i*nx+j;
      double dx = xn[nIndex] - xn[nIndex+1];
      double dy = yn[nIndex+j] - yn[nIndex+1];

      int fIndex = i*C.nx_c+j;
      Ai[fIndex] = sqrt(dx*dx + dy*dy);    
      nix[fIndex] = -1.0*dy / Ai[fIndex];
      niy[fIndex] = dx / Ai[fIndex];
    }
  }
}

void extrapCopyCoords(
    double xl_g[],        // Output: x left ghost cell coords
    double xr_g[],        // Output: x right ghost cell coords
    double xb_g[],        // Output: x bottom ghost cell coords
    double xt_g[],        // Output: x top ghost cell coords
    double yl_g[],        // Output: y left ghost cell coords
    double yr_g[],        // Output: y right ghost cell coords
    double yb_g[],        // Output: y bottom ghost cell coords
    double yt_g[],        // Output: y top ghost cell coords
    double xc[],          // Input: x center coords
    double yc[],          // Input: y center coords
    constants C)          // Input: constants including num_ghost
{
  if(C.num_ghost < 2)
  {
    cerr << "ERROR: Need to use 2 or more Ghost Layers" << endl;
    exit(-1);
  }
  for( int i=0; i < C.ny_c; i++)
  {
    // Left 
    xl_g[i*C.num_ghost] = 2.0*xc[i*C.nx_c] - xc[i*C.nx_c+1];
    xl_g[i*C.num_ghost+1] = 2.0*xl_g[i*C.num_ghost] - xc[i*C.nx_c];
    yl_g[i*C.num_ghost] = 2.0*yc[i*C.nx_c] - yc[i*C.nx_c+1];
    yl_g[i*C.num_ghost+1] = 2.0*yl_g[i*C.num_ghost] - yc[i*C.nx_c];

    // Right
    xr_g[i*C.num_ghost] = 2.0*xc[(i+1)*(C.nx_c)-1] - xc[(i+1)*(C.nx_c)-2];
    xr_g[i*C.num_ghost+1] = 2.0*xr_g[i*C.num_ghost] - xc[(i+1)*(C.nx_c)-1];
    yr_g[i*C.num_ghost] = 2.0*yc[(i+1)*(C.nx_c)-1] - yc[(i+1)*(C.nx_c)-2];
    yr_g[i*C.num_ghost+1] = 2.0*yr_g[i*C.num_ghost] - yc[(i+1)*(C.nx_c)-1];

    // Extrapolate the rest 
    for( int j=2; j < C.num_ghost; j++)
    {
      xr_g[i*C.num_ghost+j] = 2.0*xr_g[i*C.num_ghost+j-1] -xr_g[i*C.num_ghost+j-2];
      yr_g[i*C.num_ghost+j] = 2.0*yr_g[i*C.num_ghost+j-1] -yr_g[i*C.num_ghost+j-2];
      xl_g[i*C.num_ghost+j] = 2.0*xl_g[i*C.num_ghost+j-1] -xl_g[i*C.num_ghost+j-2];
      yl_g[i*C.num_ghost+j] = 2.0*yl_g[i*C.num_ghost+j-1] -yl_g[i*C.num_ghost+j-2];
    }
  }
  // top bottom
  for(int j=0; j < C.nx_c; j++)
  {
    // top 
    xt_g[j*C.num_ghost] = 2.0*xc[C.nx_c*(C.ny_c-1)+j] - xc[C.nx_c*(C.ny_c-2)+j];
    xt_g[j*C.num_ghost+1] = 2.0*xt_g[j*C.num_ghost] - xc[C.nx_c*(C.ny_c-1)+j] ;
    yt_g[j*C.num_ghost] = 2.0*yc[C.nx_c*(C.ny_c-1)+j] - yc[C.nx_c*(C.ny_c-2)+j];
    yt_g[j*C.num_ghost+1] = 2.0*yt_g[j*C.num_ghost] - yc[C.nx_c*(C.ny_c-1)+j] ;

    // bottom 
    xb_g[j*C.num_ghost] = 2.0*xc[j] - xc[C.nx_c+j];
    xb_g[j*C.num_ghost+1] = 2.0*xb_g[j*C.num_ghost] - xc[j] ;
    yb_g[j*C.num_ghost] = 2.0*yc[j] - yc[C.nx_c+j];
    yb_g[j*C.num_ghost+1] = 2.0*yb_g[j*C.num_ghost] - yc[j] ;

    for(int i=2; i < C.num_ghost; i++)
    {
      xt_g[j*C.num_ghost+i] = 2.0*xt_g[j*C.num_ghost+i-1] -xt_g[j*C.num_ghost+i-2];
      yt_g[j*C.num_ghost+i] = 2.0*yt_g[j*C.num_ghost+i-1] -yt_g[j*C.num_ghost+i-2];
      xb_g[j*C.num_ghost+i] = 2.0*xb_g[j*C.num_ghost+i-1] -xb_g[j*C.num_ghost+i-2];
      yb_g[j*C.num_ghost+i] = 2.0*yb_g[j*C.num_ghost+i-1] -yb_g[j*C.num_ghost+i-2];
    }
  }
}


void computeRTheta(
    double rc[],          // Output: Radius coord at center coord
    double thetc[],       // Output: Theta coord
    double xc[],          // Output: x center coord 
    double yc[],          // Output: x center coord 
    constants C)          // Input: number of center in x 
{
  int index; double x, y;
  for(int i = 0; i < C.ny_c; i++){
    for(int j = 0; j < C.nx_c; j++)
    {
      index = i*C.nx_c+j;
      x = xc[index]; y = yc[index];
      rc[index] = sqrt(x*x+y*y);
      thetc[index] = atan(y/x);
    }
  }
}

void getCoord(
    double xn[],          // Output: x node coord
    double yn[],          // Output: y node coord
    double xc[],          // Output: x center coord
    double yc[],          // Output: y center coord
    constants C)          // Input: Number of cells in dirs
{
  int nx = C.nx_c+1;
  int ny = C.ny_c+1;

  string filename = readMeshName(C);
  const char *FILENAME = filename.c_str();
  //open the instream
  ifstream infile;
  infile.open(FILENAME);

  double readval;
  vector<double> data;
  string line;

  if (infile.is_open())
  {
    while ( getline (infile,line))
    {
      istringstream streamA(line);
      while (streamA >> readval)
      {
        data.push_back(readval);
      }
    }
    infile.close();
  }

  for (int i=0; i < nx*ny; i++)
  {
    int datindex = i+3;
    xn[i] = data[datindex];
    yn[i] = data[datindex + nx*ny];
  }

  for (int i=0; i < C.ny_c; i++)
  {
    for (int j=0; j < C.nx_c; j++)
    {

      xc[i*C.nx_c+j] = 0.25*( xn[nx*i+j] + xn[nx*i+j+1] + xn[nx*(i+1)+j] + xn[nx*(i+1)+j+1] );
      yc[i*C.nx_c+j] = 0.25*( yn[nx*i+j] + yn[nx*i+j+1] + yn[nx*(i+1)+j] + yn[nx*(i+1)+j+1] );
    }
  }
}

void meshSize(
    int* nx_i,            // Output: Number of nodes x dir
    int* ny_i,            // Output: Number of nodes y dir
    constants C)          // Input: Constants points to which mesh
{
  string filename = readMeshName(C);
  const char *FILENAME = filename.c_str();
  //open the instream
  ifstream infile;
  infile.open(FILENAME);

  int readval;
  vector<int> data;
  string line;
  if (infile.is_open())
  {
    int i=0;
    while ( getline (infile,line) && i < 3)
    {
      istringstream streamA(line);
      while (streamA >> readval)
      {
        data.push_back(readval);
      }
      i++;
    }
    infile.close();
  }

  *nx_i = data[1];
  *ny_i = data[2];

}

string readMeshName(
    // This function reads in the mesh name and return the filename
    // it will be used to extract the data
    constants C
    ) 
{   

  if (C.f_mesh == 1) // read in smallest mesh
    return "mesh/debugMatlab.msh";
  if (C.f_mesh == 2) 
    return "grids/Curvilinear/2d17.grd";

  return "failure"; // return failure because it can not recognize the file
}

void outputArray(
    string Address,       // Output: Folder location
    string FileName,      // Output: File name
    double out[],         // Input: array to output double only
    int size,             // Input: size of array because its not passed
    int n)                // Input: Iteration number
{
  ostringstream StrConvert;
  StrConvert << n; // import n as an ostringstream
  string num_iter = StrConvert.str();// convert to string
  string Suffix = ".txt";
  string Sub = "-";
  Address = Address + "/" + FileName + Sub + num_iter + Suffix; // can combine string
  ofstream outfile; // output
  outfile.open(Address.c_str()); // access string component

  for(int i = 0; i < size; i++)
    outfile << setprecision(14) << out[i]  << endl;

  outfile.close();
}

double SIGN(double a, double b)
{
  if (b<0)
    return -a;
  else
    return a;
}

/*

//printf("ghostInex = %d, bInt = %d,  tInt= %d, bFnt = %d, tFnt = %d\n",g, bInt, tInt, bFnt, tFnt) ;
//printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xr_g[i*C.num_ghost+1],xr_g[i*C.num_ghost] , xc[(i+1)*(C.nx_c)-1], xc[(i+1)*(C.nx_c)-2] );

//printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",yb_g[j*C.num_ghost+1],yb_g[j*C.num_ghost] , yc[j], yc[C.nx_c+j] );
//printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xb_g[j*C.num_ghost+1],xb_g[j*C.num_ghost] , xc[C.nx_c*(C.ny_c-1)+j], xc[C.nx_c*(C.ny_c-2)+j] );
//printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n", j*C.num_ghost,C.nx_c*(C.ny_c-1)+j ,C.nx_c*(C.ny_c-2)+j) ;


cout << nx*i+j <<" " << nx*i+j+1  <<  " " << nx*(i+1)+j  << " " <<nx*(i+1)+j+1  << endl;
cout << xc[i*C.nx_c+j] <<" " << yc[i*C.nx_c+j] <<  " " << xn[nx*i+j] << endl;
cout << C.nx_c << " " << C.ny_c<< endl;
printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xr_g[i*C.num_ghost+1],xr_g[i*C.num_ghost] , xc[(i+1)*(C.nx_c)-1], xc[(i+1)*(C.nx_c)-2] );
printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n",i*C.num_ghost ,  (i+1)*C.nx_c-1, (i+1)*C.nx_c-2);
printf("ghostInex = %d, xcInt1st = %d,  xcInt2nd = %d,\n",i*C.num_ghost ,  i*C.nx_c, i*C.nx_c+1);
printf("ghost2 = %lf, ghost1 = %lf, xcInt1st = %lf,  xcInt2nd = %lf,\n",xl_g[i*C.num_ghost+1],xl_g[i*C.num_ghost] ,  xc[i*C.nx_c], xc[i*C.nx_c+1]);
printf("xc1stInt = %lf, and xc2ndInt = %lf \n", xc[i*C.ny_c], xc[i*C.ny_c+1]);
xl_g[i*C.num_ghost] = 2.0*xc[i*C.nx_c+1] - xc[i*C.nx_c+2];
xl_g[i*C.num_ghost+1] = 2.0*xl_g[i*C.num_ghost] - xc[i*C.nx_c+1];

printf("xl1st = %lf, and xl2nd = %lf \n", xl_g[i*C.num_ghost], xl_g[i*C.num_ghost+1]);

//cout << Aj[i*nx+j] << endl;
//cout << dx <<" "  << dy  <<  endl;
//printf("faceIndex = %d, BotNode = %d,  TopNode = %d,\n", i*nx+j ,i*nx+j , (i+1)*nx+j) ;

//cout << dx << " " << dy << endl;
//cout << Ai[i*C.nx_c+j] << endl;
//printf("faceIndex = %d, leftNode = %d,  rightNode = %d,\n", i*C.nx_c+j, i*nx+j, i*nx+j+1 ) ;

*/


#include "mhdRT.hpp"


// constructor


Map2Eigen :: Map2Eigen(int nj, int ni, int nequ):
    Q_raw(new double[nj*ni*nequ]),
      Q{{Q_raw,   nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+1, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+2, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+3, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+4, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+5, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+6, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) },
        {Q_raw+7, nj, ni, Stride<Dynamic,Dynamic>(ni*nequ,nequ) } }
{
  
}

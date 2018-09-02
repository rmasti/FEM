#include "mhdRT.hpp"


// constructor


Map2Eigen :: Map2Eigen(int ni, int nj, int nequ):
    Q_raw(new double[ni*nj*nequ]),
      Q{ {Q_raw, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+1, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+2, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+3, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+4, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+5, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+6, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) },
        {Q_raw+7, ni, nj, Stride<Dynamic,Dynamic>(ni*nequ, nequ) } }
{
  
}

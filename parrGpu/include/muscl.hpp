double limiter(double r, const int* lim);

double limiter(const double r, const int& lim)
{
  double theta;
  switch (lim)
  {
    case 1:
      theta = 1.0;
    case 2:
      theta = (r+fabs(r))/(1.0+fabs(r));
    case 3:
      theta = (r+r*r)/(1.0+r*r);
    case 4:
      theta = 1.5*(r*r + r) / (1.0 + r + r*r);
    case 5:
      theta = fmax(0,fmin(2.0*r,fmin(0.5*(1.0+r),2.0)));
    case 6:
      theta = fmax(0,fmin(1.0,r));
    case 7:
      theta = fmax(0,fmax(fmin(2.0*r,1.0),fmin(r,2.0)));
  }
  return theta;
}

#include <math.h>
#include "response.h"

// ===============================
// Function: Hyperbolic Arctangent
// ===============================
double atanh(double x)
{
  double r;

  r=.5*log((1.0+x)/(1.0-x));

  return r;
}

// =======================
// Function: concentration
// =======================
double conc(double a, double b, double p,double q)
{
  double
    r;

  r = (1.0+a*p)/(1.0-a*p)+(1.0-b*p)*atanh(2*q-1.0)/(1.0+b*p);

  return r;
}
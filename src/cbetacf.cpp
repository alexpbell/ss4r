/**
@file
@author David Fournier
@copyright Copyright (c) 2009-2020 ADMB Foundation

@brief Wrapper for betacf function.
*/
#include "fvar.h"
#include "betacf_val.h"

double betacf(const double a, const double b, const double x, int MAXIT){
  typedef double Float;
  return betacf<Float>(a,b,x,MAXIT);
}

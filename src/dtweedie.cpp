#include <vector>
#include <type_traits>
#include "fvar.h"

#ifdef DEBUG
  #include <cassert>
  #include <climits>
#endif

unsigned int imax2(unsigned int a, double v)
{
#ifdef DEBUG
  assert(v <= double(UINT_MAX));
#endif
  return static_cast<double>(a) >= v ? a : static_cast<unsigned int>(v);
}
unsigned int imin2(unsigned int a, unsigned int b)
{
  return a < b ? a : b;  
}
double fmax2(double a, const dvariable& v)
{
  double b = value(v);
  return a > b ? a : b;  
}
double asDouble(const dvariable& v)
{
  return value(v);
}
dvariable lgamma(const prevariable& v)
{
  return gammln(v);  
}

/**
  \brief dtweedie function (same as dtweedie.series from R package
  'tweedie').

  Silently returns NaN if not within the valid parameter range:
  \f[ (0 \leq y) \land (0 < \mu) \land (0 < \phi) \land (1 < p) \land (p < 2) \f] .

  \note Parameter order differs from the R version.

  \warning The derivative wrt. the y argument is disabled
  (zero). Hence the tweedie distribution can only be used for *data*
  (not random effects).

  The dtweedie function was copied from TMB repository at
  https://github.com/kaskr/adcomp, then modified for ADMB
  data types.

  \ingroup R_style_distribution
*/
template<class Type>
Type _dtweedie(const double y, Type& mu, Type& phi, Type& p, bool give_log) {
  Type p1 = p - 1.0, p2 = 2.0 - p;
  Type ans = -pow(mu, p2) / (phi * p2); // log(prob(y=0))
  if (y > 0) {
    ans += tweedie_logW(y, phi, p);
    ans += -y / (phi * p1 * pow(mu, p1)) - log(y);
  }
  return ( give_log ? ans : exp(ans) );
}

/// dtweedie is a wrapper to _dtweedie using ADMB data types. 
dvariable dtweedie(const double y, dvariable& mu, dvariable& phi, dvariable& p, const bool use_log)
{
  return _dtweedie(y, mu, phi, p, use_log);
}
// Re-structured version of tweedie density function from 'cplm' package.

/* the threshold used in finding the bounds of the series */
#define TWEEDIE_DROP 37.0
/* the loop increment used in finding the bounds of the series */
#define TWEEDIE_INCRE 5
/* the max number of terms allowed in the finite sum approximation*/
#define TWEEDIE_NTERM 20000

/** \brief Calculate \f$\log W(y, \phi, p)$\f with notation as in Dunn
    and Smyth 2005 page 269 equation 2.  Required to calculate the
    density of the Tweedie distribution.

    The tweedie_logW function was copied from TMB repository at 
    https://github.com/kaskr/adcomp, then modified for ADMB
    data types.

    \param y _positive_ observation
    \param phi scalar: the dispersion parameter
    \param p scalar: the index parameter
*/
template<class Float>
Float tweedie_logW(double y, Float& phi, Float& p){
  bool ok = (0 < y) && (0 < phi) && (1 < p) && (p < 2);
  if (!ok) return NAN;

  Float p1 = p - 1.0, p2 = 2.0 - p;
  Float a = - p2 / p1, a1 = 1.0 / p1;
  Float cc, w, sum_ww = 0.0;
  double ww_max = -INFINITY ;

  /* only need the lower bound and the # terms to be stored */
  double jmax = 0;
  Float logz = 0;

  /* compute jmax for the given y > 0*/
  cc = a * log(p1) - log(p2);
  jmax = asDouble( fmax2(1.0, pow(y, p2) / (phi * p2)) );
  logz = - a * log(y) - a1 * log(phi) + cc;

  /* find bounds in the summation */
  /* locate upper bound */
  cc = logz + a1 + a * log(-a);
  double j = jmax ;
  w = a1 * j ;
  while (1) {
    j += TWEEDIE_INCRE ;
    if (j * (cc - a1 * log(j)) < (w - TWEEDIE_DROP))
      break ;
  }
  unsigned int jh = static_cast<unsigned int>(ceil(j));
  /* locate lower bound */
  j = jmax;
  while (1) {
    j -= TWEEDIE_INCRE ;
    if (j < 1 || j * (cc - a1 * log(j)) < w - TWEEDIE_DROP)
      break ;
  }
  unsigned int jl = imax2(1, floor(j)) ;
  unsigned int jd = jh - jl + 1;

  /* set limit for # terms in the sum */
  unsigned int nterms = imin2(jd, TWEEDIE_NTERM) ;
  //std::vector<Float> ww(nterms);
  /* evaluate series using the finite sum*/
  /* y > 0 */
  sum_ww = 0.0 ;
  unsigned int iterm = imin2(jd, nterms) ;
  std::vector<Float> ww(iterm);
  for (unsigned int k = 0; k < iterm; k++) {
    j = k + jl ;
    ww[k] = j * logz - lgamma(1 + j) - lgamma(-a * j);
    ww_max = fmax2( ww_max, asDouble(ww[k]) );
  }
  for (unsigned int k = 0; k < iterm; k++)
    sum_ww += exp(ww[k] - ww_max);
  Float ans = log(sum_ww) + ww_max  ;

  return ans;
}

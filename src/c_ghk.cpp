/**
@file
@author David Fournier
@copyright Copyright (c) 2008-2020 Regents of the University of California
*/
#include "fvar.h"

/**
 * Description not yet available.
 * \param
 */
double ghk(const dvector& lower,const dvector& upper,const dmatrix& Sigma,
  const dmatrix& eps)
{
  int m=eps.indexmax();
  int n=lower.indexmax();
  double ssum=0.0;
  dmatrix ch=choleski_decomp(Sigma);
  dvector l(1,n);
  dvector u(1,n);

  for (int k=1;k<=m;k++)
  {
    double weight=1.0;
    l=lower;
    u=upper;
    for (int j=1;j<=n;j++)
    {
      l(j)/=ch(j,j);
      u(j)/=ch(j,j);
      double Phiu=cumd_norm(u(j));
      double Phil=cumd_norm(l(j));
      weight*=Phiu-Phil;
      double eta=inv_cumd_norm((Phiu-Phil)*eps(k,j)+Phil);
      for (int i=j+1;i<=n;i++)
      {
        double tmp=ch(i,j)*eta;
        l(i)-=tmp;
        u(i)-=tmp;
      }
    }
    ssum+=weight;
  }
  return ssum/m;
}

/**
 * Description not yet available.
 * \param
 */
void ghk_test(const dmatrix& eps,int i)
{
  if (i<eps.indexmin())
  {
    ad_exit(21);
  }
  else  if (i>eps.indexmax())
  {
    ad_exit(21);
  }
}

/**
 * Description not yet available.
 * \param
 */
double ghk(const dvector& lower,const dvector& upper,const dmatrix& Sigma,
  const dmatrix& eps,int _i)
{
  int n=lower.indexmax();
  dmatrix ch=choleski_decomp(Sigma);
  dvector l(1,n);
  dvector u(1,n);

  ghk_test(eps,_i); // test for valid i range
  double weight=1.0;
  int k=_i;
  {
    l=lower;
    u=upper;
    for (int j=1;j<=n;j++)
    {
      l(j)/=ch(j,j);
      u(j)/=ch(j,j);
      double Phiu=cumd_norm(u(j));
      double Phil=cumd_norm(l(j));
      weight*=Phiu-Phil;
      double eta=inv_cumd_norm((Phiu-Phil)*eps(k,j)+Phil);
      for (int i=j+1;i<=n;i++)
      {
        double tmp=ch(i,j)*eta;
        l(i)-=tmp;
        u(i)-=tmp;
      }
    }
  }
  return weight;
}

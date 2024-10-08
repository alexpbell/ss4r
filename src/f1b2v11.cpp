/*
Author: David Fournier
Copyright (c) 2008-2012 Regents of the University of California
*/

#include "df1b2fun.h"

/**
Returns a df1b2vector containing the differences of an _x(i) and  _x(i + 1) for i = 1 to x.indexmax() - 1.

\param _x input.
*/
df1b2vector first_difference(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  df1b2vector tmp;
  tmp.noallocate(mmin,mmax-1);
  int i;
  for (i=mmin;i<mmax;i++)
  {
    tmp(i)=x(i+1)-x(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector square(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  df1b2vector tmp;
  tmp.noallocate(mmin,mmax);
  int i;
  for (i=mmin;i<=mmax;i++)
  {
    tmp(i)=square(x(i));
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector cube(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  df1b2vector tmp;
  tmp.noallocate(mmin,mmax);
  int i;
  for (i=mmin;i<=mmax;i++)
  {
    tmp(i)=square(x(i));
  }
  return tmp;
}

/*
 * $Id$
 *
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
/**
 * \file
 * Description not yet available.
 */
#include "fvar.h"

/**
 * Description not yet available.
 * \param
 */
dvector dvector::operator-(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector tmp(mmin,mmax);
  double* pvi = v + mmin;
  double* ptmpi = tmp.get_v() + mmin;
  for (int i=mmin;i<=mmax;i++)
  {
    *ptmpi = -(*pvi);

    ++pvi;
    ++ptmpi;
  }
  return tmp;
}

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
// file fvar.cpp
// constructors, destructors and misc functions involving class dvariable

#include "fvar.h"

#ifdef __TURBOC__
  #pragma hdrstop
  #include <iostream.h>
#endif

#ifdef __ZTC__
  #include <iostream.hpp>
#endif


#include <stdio.h>
#ifndef __SUN__
#endif
#include <math.h>

/**
 * Description not yet available.
 * \param
 */
dmatrix operator/(const dmatrix& m, const double e)
{
  dmatrix tmp;
  tmp.allocate(m);
  for (int i=m.rowmin();i<=m.rowmax();i++)
  {
    tmp(i)=m(i)/e;
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
dmatrix operator/(const double e, const dmatrix& m)
{
  dmatrix tmp;
  tmp.allocate(m);
  for (int i=m.rowmin();i<=m.rowmax();i++)
  {
    tmp(i)=e/m(i);
  }
  return tmp;
}

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
// file: dmat_io.cpp

// i/o operations for class dmatrix
#include "fvar.h"

#ifdef __TURBOC__
  #pragma hdrstop
  #include <iostream.h>
  #include <iomanip.h>
  #include <fstream.h>
  #define __USE_IOSTREAM__
#endif

#ifdef __ZTC__
  #include <iostream.hpp>
  #include <iomanip.hpp>
  #include <fstream.hpp>
  #define __USE_IOSTREAM__
#endif

#include <string.h>

/**
 * Description not yet available.
 * \param
 */
ostream& operator<<(const ostream& _ostr, const imatrix& z)
{
  ostream& ostr= (ostream&) _ostr;
  int min=z.rowmin();
  int max=z.rowmax();
  for (int i=min;i<=max-1;i++)
  {
     ostr << z[i] << endl;
  }
  ostr << z[max];
  return ostr;
}

/**
 * Description not yet available.
 * \param
 */
istream& operator>>(const istream& _istr,const imatrix& z)
{
  istream& istr= (istream&) _istr;
  int min=z.rowmin();
  int max=z.rowmax();
  for (int i=min;i<=max;i++)
  {
     istr >> z[i];
  }
  return istr;
}

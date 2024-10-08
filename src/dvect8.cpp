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

#ifdef DEBUG
  #include <cassert>
  #include <climits>
#endif

/**
 * Description not yet available.
 * \param
 */
dvector::dvector(const ivector& u)
 {
   allocate(u.indexmin(),u.indexmax());
   for ( int i=indexmin(); i<=indexmax(); i++)
   {
     elem(i)=u.elem(i);
   }
 }

/**
Construct dvector with values in lvector.

\param u values to copy
*/
dvector::dvector(const lvector& u)
{
  allocate(u.indexmin(), u.indexmax());
  for (int i=indexmin(); i <= indexmax(); i++)
  {
    elem(i) = static_cast<double>(u.elem(i));
  }
}

/**
 * Description not yet available.
 * \param
 */
dvector dvector::operator ()(const ivector& u)
 {
   dvector tmp(u.indexmin(),u.indexmax());

   for ( int i=u.indexmin(); i<=u.indexmax(); i++)
   {
     tmp(i)=(*this)(u(i));
   }
   return tmp;
 }

/**
 * Description not yet available.
 * \param
 */
dvector dvector::operator ()(const lvector& u)
 {
   dvector tmp(u.indexmin(),u.indexmax());

   for ( int i=u.indexmin(); i<=u.indexmax(); i++)
   {
#ifdef DEBUG
     assert(u(i) <= INT_MAX);
#endif
     tmp(i)=(*this)(static_cast<int>(u(i)));
   }
   return tmp;
 }

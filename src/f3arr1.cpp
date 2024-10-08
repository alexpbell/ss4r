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
#include "admb_messages.h"

#ifndef OPT_LIB

/**
 * Description not yet available.
 * \param
 */
prevariable dvar3_array::operator () (int k, int i, int j)
{
  if (!allocated(*this))
  {
    cerr << "trying to access an unallocated object" << endl;
    ad_exit(21);
  }
  if (k < slicemin())
  {
    ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too low",
      "prevariable dvar3_array::operator () (int k, int i, int j)",
    slicemin(), slicemax(), k);
  }
  if (k > slicemax())
  {
    ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too high",
      "prevariable dvar3_array::operator () (int k, int i, int j)",
    slicemin(), slicemax(), k);
  }
  return (t + k)->operator()(i,j);
}

/**
 * Description not yet available.
 * \param
 */
dvar_vector& dvar3_array::operator () (int k, int i)
{
  if (k < slicemin())
  {
    ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too low",
      "dvar_vector& dvar3_array::operator () (int k, int i)",
    slicemin(), slicemax(), k);
  }
  if (k > slicemax())
  {
    ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too high",
      "dvar_vector& dvar3_array::operator () (int k, int i)",
    slicemin(), slicemax(), k);
  }
  return (t + k)->operator()(i);
}

/**
 * Description not yet available.
 * \param
 */
 dvar_matrix& dvar3_array::operator[] (int i)
 {
   if (i < slicemin())
   {
     ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too low",
     "dvar_matrix& dvar3_array::operator [] (int i)",
     slicemin(), slicemax(), i);
   }
   if (i > slicemax())
   {
     ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too high",
     "dvar_matrix& dvar3_array::operator [] (int i)",
     slicemin(), slicemax(), i);
   }
   return t[i];
 }

/**
 * Description not yet available.
 * \param
 */
 dvar_matrix& dvar3_array::operator() (int i)
 {
   if (!allocated(*this))
   {
       cerr << "trying to access an unallocated object" << endl;
       ad_exit(21);
   }
   if (i < slicemin())
   {
     ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too low",
     "dvar_matrix& dvar3_array::operator () (int i)",
     slicemin(), slicemax(), i);
   }
   if (i > slicemax())
   {
     ADMB_ARRAY_BOUNDS_ERROR("array bound exceeded -- slice index too high",
     "dvar_matrix& dvar3_array::operator () (int i)",
     slicemin(), slicemax(), i);
   }
   return t[i];
 }
#endif

/**
 * Description not yet available.
 * \param
 */
dvariable sum(const dvar3_array& m)
{
  gradient_structure* gs = gradient_structure::_instance;
  gs->RETURN_ARRAYS_INCREMENT();

  int min = m.indexmin();
  int max = m.indexmax();
  dvariable tmp=0.;
  const dvar_matrix* pmi = &m(min);
  for (int i = min; i <= max; ++i)
  {
    tmp += sum(*pmi);
    ++pmi;
  }
  gs->RETURN_ARRAYS_DECREMENT();
  return tmp;
}

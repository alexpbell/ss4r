/*
 * $Id$
 *
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "admodel.h"

 void param_init_number_vector::set_initial_value(const double_index_type& _it)
 {
    it=new double_index_type(_it);
 }
/**
Default constructor
*/
param_init_number_vector::param_init_number_vector():
  v(NULL),
  index_min(0),
  index_max(0),
  it(NULL)
{
}
/**
Destructor
*/
param_init_number_vector::~param_init_number_vector()
{
  deallocate();
}
/**
Free member allocated memory.
*/
void param_init_number_vector::deallocate(void)
{
  if(it)
  {
    delete it;
    it = NULL;
  }
  if (v)
  {
    v+=indexmin();
    delete [] v;
    v = NULL;
  }
}

void param_init_number_vector::allocate(
  int min1,
  int max1,
  const char* s)
{
  allocate(min1,max1,1,s);
}

void param_init_number_vector::allocate(
  int min1,
  int max1,
  const index_type& phase_start,
  const char* s)
{
  int size = max1 - min1 + 1;
  if (size > 0)
  {
    v = new param_init_number[static_cast<unsigned int>(size)];
    if (!v)
    {
        cerr << " error trying to allocate memory in "
          "param_init_vector_vector " << endl;
        ad_exit(1);
    }

    index_min=min1;
    index_max=max1;
    v-=indexmin();
    for (int i=indexmin();i<=indexmax();i++)
    {
       if (it) v[i].set_initial_value(ad_double((*it)[i]));
       adstring ss=s + adstring("[") + str(i) + adstring("]");
       v[i].allocate(ad_integer(phase_start[i]),(char*)(ss) );
    }
  }
}

  dvector value(const param_init_number_vector& _t)
  {
    ADUNCONST(param_init_number_vector, t);
    const int min = t.indexmin();
    const int max = t.indexmax();
    dvector vt(min, max);
    for (int i = min; i <= max; i++)
    {
       vt(i) = value(t(i));
    }
    return vt;
  }

  dvector value(const param_init_bounded_number_vector& _t)
  {
    ADUNCONST(param_init_bounded_number_vector, t);
    const int min = t.indexmin();
    const int max = t.indexmax();
    dvector vt(min, max);
    for (int i = min; i <= max; i++)
    {
       vt(i) = value(t(i));
    }
    return vt;
  }

/**
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "admodel.h"

/**
Sets the initial value for elements of param_init_bounded_vector.

\param initial_value double
*/
void param_init_bounded_number_vector::set_initial_value(
  const double_index_type& initial_value)
{
  if (it)
  {
    delete it;
    it = 0;
  }
  it = new double_index_type(initial_value);

  if (v)
  {
    int mmin = indexmin();
    int mmax = indexmax();
    for (int i = mmin; i <= mmax; ++i)
    {
      v[i] = ad_double(initial_value(i));
    }
  }
}
/// Default constructor
param_init_bounded_number_vector::param_init_bounded_number_vector():
  v(nullptr),
  index_min(1),
  index_max(0),
  it(nullptr)
{
}
/// Destructor
param_init_bounded_number_vector::~param_init_bounded_number_vector()
{
  deallocate();
}
/// Deallocate param_init_bounded_number_vector memory.
void param_init_bounded_number_vector::deallocate(void)
{
  if (it)
  {
    delete it;
    it = NULL;
  }
  if (v)
  {
    v += indexmin();
    delete [] v;
    v = NULL;
  }
  index_min = 1;
  index_max = 0;
}
/**
Overload the allocate function to use a data_matrix object.
\author Steve Martell
*/
void param_init_bounded_number_vector::allocate(const data_matrix &m,
  const char *s)
{
  int min1 = m.rowmin();
  int max1 = m.rowmax();
  double_index_type bmin = column(m,1);
  double_index_type bmax = column(m,2);
  index_type phz1 = ivector(column(m,3));
  allocate(min1,max1,bmin,bmax,phz1,s);
}
/**
Allocate vector of param_init_bounded_number with dimension [min1 to max1]
and bounded by [bmin, bmax].

Note: phase_start is defaulted to 1.

\param min1 lower vector index
\param max1 upper vector index
\param bmin lower bounds
\param bmax upper bounds
\param s id
*/
void param_init_bounded_number_vector::allocate(
  int min1,
  int max1,
  const double_index_type& bmin,
  const double_index_type& bmax,
  const char* s)
{
  allocate(min1, max1, bmin, bmax, 1, s);
}
/**
Allocate vector of param_init_bounded_number with dimension [min1 to max1]
and bounded by [bmin, bmax].

\param min1 lower vector index
\param max1 upper vector index
\param bmin lower bounds
\param bmax upper bounds
\param phase_start
\param s id
*/
void param_init_bounded_number_vector::allocate(
  int min1,
  int max1,
  const double_index_type& bmin,
  const double_index_type& bmax,
  const index_type& phase_start,
  const char* s)
{
  if (max1 >= min1)
  {
    v = new param_init_bounded_number[
          static_cast<unsigned int>(max1 - min1 + 1)];
    if (!v)
    {
      cerr << " error trying to allocate memory in "
           << "param_init_bounded_number_vector " << endl;
      ad_exit(1);
    }
    index_min = min1;
    index_max = max1;
    v -= indexmin();
    for (int i = indexmin(); i <= indexmax(); ++i)
    {
      if (it) v[i].set_initial_value(ad_double((*it)[i]));
      adstring ss=s + adstring("[") + str(i) + adstring("]");
      v[i].allocate(ad_double(bmin[i]),ad_double(bmax[i]),
        ad_integer(phase_start[i]),(char*)(ss) );
    }
  }
}
/**
Returns vector of scalefactor values.
*/
dvector param_init_number_vector::get_scalefactor()
{
  int mmin = indexmin();
  int mmax = indexmax();
  dvector scalefactors(mmin, mmax);
  for (int i = mmin; i <= mmax; ++i)
  {
    scalefactors(i) = (*this)(i).get_scalefactor();
  }
  return scalefactors;
}
/**
Sets individual scalefactor for param_init_number_vector.

\param scalefactor dvector
*/
void param_init_number_vector::set_scalefactor(const dvector& scalefactors)
{
  int mmin = indexmin();
  int mmax = indexmax();
  if (scalefactors.indexmin() != mmin || scalefactors.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
         << " init_number_vector::set_scalefactor" << endl;
    ad_exit(1);
  }
  for (int i = mmin; i <= mmax; ++i)
  {
    (*this)(i).set_scalefactor(scalefactors(i));
  }
}
/**
Sets scalefactor for param_init_number_vector.

\param scalefactor double
*/
void param_init_number_vector::set_scalefactor(double scalefactor)
{
  int mmin = indexmin();
  int mmax = indexmax();
  for (int i = mmin; i <= mmax; ++i)
  {
    (*this)(i).set_scalefactor(scalefactor);
  }
}

dvector param_init_vector_vector::get_scalefactor(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector s(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    s(i)=(*this)(i).get_scalefactor();
  }
  return s;
}
void param_init_vector_vector::set_scalefactor(const dvector& s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  if (s.indexmin()!=mmin || s.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
     " init_vector_vector::set_scalefactor" << endl;
    ad_exit(1);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s(i));
  }
}

void param_init_vector_vector::set_scalefactor(double s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s);
  }
}

dvector param_init_matrix_vector::get_scalefactor(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector s(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    s(i)=(*this)(i).get_scalefactor();
  }
  return s;
}
void param_init_matrix_vector::set_scalefactor(const dvector& s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  if (s.indexmin()!=mmin || s.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
     " init_matrix_vector::set_scalefactor" << endl;
    ad_exit(1);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s(i));
  }
}

void param_init_matrix_vector::set_scalefactor(double s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s);
  }
}

dvector param_init_bounded_number_vector::get_scalefactor(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector s(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    s(i)=(*this)(i).get_scalefactor();
  }
  return s;
}
void param_init_bounded_number_vector::set_scalefactor(const dvector& s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  if (s.indexmin()!=mmin || s.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
     " init_bounded_number_vector::set_scalefactor" << endl;
    ad_exit(1);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s(i));
  }
}

void param_init_bounded_number_vector::set_scalefactor(double s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s);
  }
}

dvector param_init_bounded_vector_vector::get_scalefactor(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector s(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    s(i)=(*this)(i).get_scalefactor();
  }
  return s;
}
void param_init_bounded_vector_vector::set_scalefactor(const dvector& s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  if (s.indexmin()!=mmin || s.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
     " init_bounded_vector_vector::set_scalefactor" << endl;
    ad_exit(1);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s(i));
  }
}

void param_init_bounded_vector_vector::set_scalefactor(double s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s);
  }
}

dvector param_init_bounded_matrix_vector::get_scalefactor(void)
{
  int mmin=indexmin();
  int mmax=indexmax();
  dvector s(mmin,mmax);
  for (int i=mmin;i<=mmax;i++)
  {
    s(i)=(*this)(i).get_scalefactor();
  }
  return s;
}
void param_init_bounded_matrix_vector::set_scalefactor(const dvector& s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  if (s.indexmin()!=mmin || s.indexmax() != mmax)
  {
    cerr << "non matching vector bounds in"
     " init_bounded_matrix_vector::set_scalefactor" << endl;
    ad_exit(1);
  }

  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s(i));
  }
}

void param_init_bounded_matrix_vector::set_scalefactor(double s)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    (*this)(i).set_scalefactor(s);
  }
}

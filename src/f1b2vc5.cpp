/**
\file
Author: David Fournier
Copyright (c) 2008-2012 Regents of the University of California
*/
#include "df1b2fun.h"
#ifdef DEBUG
  #include <cassert>
#endif

/**
 * Description not yet available.
 * \param
 */
void check_shape(const df1b2vector & _x,const dvector & _y,const char * s)
{
  ADUNCONST(df1b2vector,x)
  ADUNCONST(dvector,y)
  if (x.indexmin() != y.indexmin() || x.indexmax() != y.indexmax())
  {
    cerr << "Incompatible shapes in df1b2vector function" << s << endl;
    ad_exit(1);
  }
}

/**
 * Description not yet available.
 * \param
 */
void check_shape(const df1b2vector & _x,const df1b2vector & _y,const char * s)
{
  ADUNCONST(df1b2vector,x)
  ADUNCONST(df1b2vector,y)
  if (x.indexmin() != y.indexmin() || x.indexmax() != y.indexmax())
  {
    cerr << "Incompatible shapes in df1b2vector function" << s << endl;
    ad_exit(1);
  }
}

/**
 * Description not yet available.
 * \param
 */
void check_shape(const dvector & _x,const df1b2vector & _y,const char * s)
{
  ADUNCONST(dvector,x)
  ADUNCONST(df1b2vector,y)
  if (x.indexmin() != y.indexmin() || x.indexmax() != y.indexmax())
  {
    cerr << "Incompatible shapes in df1b2vector function" << s << endl;
    ad_exit(1);
  }
}

/**
 * Description not yet available.
 * \param
 */
void check_shape(const df1b2vector & _x,const df1b2matrix & _y,const char * s)
{
  ADUNCONST(df1b2vector,x)
  ADUNCONST(df1b2matrix,y)
  if (x.indexmin() != y(y.indexmin()).indexmin() ||
    x.indexmax() != y(y.indexmin()).indexmax())
  {
    cerr << "Incompatible shapes in df1b2vector function" << s << endl;
    ad_exit(1);
  }
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector& df1b2vector::operator=(const df1b2vector& _x)
{
  if (allocated())
  {
    ADUNCONST(df1b2vector,x)
    check_shape(*this,x,"df1b2vector& df1b2vector::operator =");
    int mmin=x.indexmin();
    int mmax=x.indexmax();
    for (int i=mmin;i<=mmax;i++)
    {
      (*this)(i)=x(i);
    }
  }
  else
  {
    copy(_x);
  }
  return *this;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector& df1b2vector::operator = (const dvector& _x)
{
  ADUNCONST(dvector,x)
  check_shape(*this,x,"df1b2vector& df1b2vector::operator =");
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++) (*this)(i)=x(i);
  return *this;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector& df1b2vector::operator = (const df1b2variable& _x)
{
  ADUNCONST(df1b2variable,x)
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++) (*this)(i)=x;
  return *this;
}
/**
Assign value x to each element of df1b2vector.

\param x assigment value
*/
df1b2vector& df1b2vector::operator=(double x)
{
  int mmin=indexmin();
  int mmax=indexmax();
  for (int i=mmin;i<=mmax;i++) (*this)(i)=x;
  return *this;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector operator * (const dmatrix& _M,const df1b2vector& _x)
{
  ADUNCONST(dmatrix,M)
  ADUNCONST(df1b2vector,x)
  //check_shape(x,M,"operator *");
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  //int mmin=x.indexmin();
  //int mmax=x.indexmax();
  df1b2vector tmp(rmin,rmax);
  tmp.initialize();
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=M(i)*x;
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector operator * (const df1b2matrix& _M,const df1b2vector& _x)
{
  ADUNCONST(df1b2matrix,M)
  ADUNCONST(df1b2vector,x)
  //check_shape(x,M,"operator *");
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  //int mmin=x.indexmin();
  //int mmax=x.indexmax();
  df1b2vector tmp(rmin,rmax);
  tmp.initialize();
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=M(i)*x;
  }
  return tmp;
}
/**
Vector and Matrix Multiplication

\param _x is a vector with dimensions (m_min, m_max).
\param _M is a matrix with dimensions (m_min, m_max, n_min, n_max).
\return vector result of _x * _M with dimensions (n_min, n_max)
*/
df1b2vector operator*(const df1b2vector& _x, const df1b2matrix& _M)
{
  ADUNCONST(df1b2matrix,M)
  ADUNCONST(df1b2vector,x)
  //check_shape(x,M,"operator *");
  int rmin = M.indexmin();
  int cmin = M(rmin).indexmin();
  int cmax = M(rmin).indexmax();
  int mmin = x.indexmin();
  int mmax = x.indexmax();
  df1b2vector tmp(cmin, cmax);
#ifndef OPT_LIB
  tmp.initialize();
#endif
  for (int i = cmin; i <= cmax; ++i)
  {
    for (int j = mmin; j <= mmax; ++j)
    {
      tmp(i) += M(j, i) * x(j);
    }
  }
  return tmp;
}
/**
Vector and Matrix Multiplication

\param _x is a vector with dimensions (m_min, m_max).
\param _M is a matrix with dimensions (m_min, m_max, n_min, n_max).
\return vector result of _x * _M with dimensions (n_min, n_max)
*/
df1b2vector operator*(const df1b2vector& _x, const dmatrix& _M)
{
  ADUNCONST(dmatrix,M)
  ADUNCONST(df1b2vector,x)
  //check_shape(x,M,"operator *");
  int rmin = M.indexmin();
  int cmin = M(rmin).indexmin();
  int cmax = M(rmin).indexmax();
  int mmin = x.indexmin();
  int mmax = x.indexmax();
  df1b2vector tmp(cmin, cmax);
#ifndef OPT_LIB
  tmp.initialize();
#endif
  for (int i = cmin; i <= cmax; ++i)
  {
    for (int j = mmin; j <= mmax; ++j)
    {
      tmp(i) += M(j, i) * x(j);
    }
  }
  return tmp;
}
/**
Vector and Matrix Multiplication

\param _x is a vector with dimensions (m_min, m_max).
\param _M is a matrix with dimensions (m_min, m_max, n_min, n_max).
\return vector result of _x * _M with dimensions (n_min, n_max)
*/
df1b2vector operator*(const dvector& _x, const df1b2matrix& _M)
{
  ADUNCONST(df1b2matrix,M)
  ADUNCONST(dvector,x)
  //check_shape(x,M,"operator *");
  int rmin = M.indexmin();
  int cmin = M(rmin).indexmin();
  int cmax = M(rmin).indexmax();
  int mmin = x.indexmin();
  int mmax = x.indexmax();
  df1b2vector tmp(cmin, cmax);
#ifndef OPT_LIB
  tmp.initialize();
#endif
  for (int i = cmin; i <= cmax; ++i)
  {
    for (int j = mmin; j <= mmax; ++j)
    {
      tmp(i) += M(j, i) * x(j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix operator * (const df1b2matrix& _MM,const df1b2matrix& _NN)
{
  df1b2matrix& M = (df1b2matrix&)_MM;
  df1b2matrix& N = (df1b2matrix&)_NN;
  //check_shape(x,M,"operator *");
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  int kmin=N.rowmin();
  int kmax=N.rowmax();
  int cmin=N(kmin).indexmin();
  int cmax=N(kmin).indexmax();
  if (M(rmin).indexmin()!=N.indexmin() || M(rmin).indexmax()!=N.indexmax())
  {
    cerr << "incompatible matrix sizes" << endl;
    ad_exit(1);
  }
  df1b2matrix tmp(rmin,rmax,cmin,cmax);
  tmp.initialize();
  for (int i=rmin;i<=rmax;i++)
  {
    for (int j=cmin;j<=cmax;j++)
    {
      for (int k=kmin;k<=kmax;k++)
      {
        tmp(i,j)+=M(i,k)*N(k,j);
      }
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector operator * (const df1b2matrix& _M,const dvector& _x)
{
  ADUNCONST(df1b2matrix,M)
  ADUNCONST(dvector,x)
  //check_shape(x,M,"operator *");
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  int cmin=x.indexmin();
  int cmax=x.indexmax();
  df1b2vector tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
    for (int j=cmin;j<=cmax;j++)
      tmp(i)+=M(i,j)*x(j);
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_prod(const df1b2matrix& _MM,const df1b2matrix& _NN)
{
  df1b2matrix& M = (df1b2matrix&)_MM;
  df1b2matrix& N = (df1b2matrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)*N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_prod(const dmatrix& _MM,const df1b2matrix& _NN)
{
  dmatrix& M = (dmatrix&)_MM;
  df1b2matrix& N = (df1b2matrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)*N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_prod(const df1b2matrix& _MM,const dmatrix& _NN)
{
  df1b2matrix& M = (df1b2matrix&)_MM;
  dmatrix& N = (dmatrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)*N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_div(const df1b2matrix& _MM,const df1b2matrix& _NN)
{
  df1b2matrix& M = (df1b2matrix&)_MM;
  df1b2matrix& N = (df1b2matrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)/N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_div(const dmatrix& _MM,const df1b2matrix& _NN)
{
  dmatrix& M = (dmatrix&)_MM;
  df1b2matrix& N = (df1b2matrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)/N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix elem_div(const df1b2matrix& _MM,const dmatrix& _NN)
{
  df1b2matrix& M = (df1b2matrix&)_MM;
  dmatrix& N = (dmatrix&)_NN;
  int rmin=M.indexmin();
  int rmax=M.indexmax();
  df1b2matrix tmp(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    int cmin=M(i).indexmin();
    int cmax=M(i).indexmax();
    tmp(i).noallocate(cmin,cmax);
    for (int j=cmin;j<=cmax;j++)
    {
      tmp(i,j)=M(i,j)/N(i,j);
    }
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_div(const dvector& _v,const df1b2vector& _w)
{
  ADUNCONST(dvector,v)
  ADUNCONST(df1b2vector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)/w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_div(const df1b2vector& _v,const df1b2vector& _w)
{
  ADUNCONST(df1b2vector,v)
  ADUNCONST(df1b2vector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)/w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_prod(const df1b2vector& _v,const df1b2vector& _w)
{
  ADUNCONST(df1b2vector,v)
  ADUNCONST(df1b2vector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)*w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_prod(const df1b2vector& _v,const dvector& _w)
{
  ADUNCONST(df1b2vector,v)
  ADUNCONST(dvector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)*w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_prod(const dvector& _v,const df1b2vector& _w)
{
  ADUNCONST(dvector,v)
  ADUNCONST(df1b2vector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)*w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2vector elem_div(const df1b2vector& _v,const dvector& _w)
{
  ADUNCONST(df1b2vector,v)
  ADUNCONST(dvector,w)
  int rmin=v.indexmin();
  int rmax=v.indexmax();
  df1b2vector tmp;
  tmp.noallocate(rmin,rmax);
  for (int i=rmin;i<=rmax;i++)
  {
    tmp(i)=v(i)/w(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix::df1b2matrix(int nrl,int nrh,int ncl,int nch)
{
  if (nrl>nrh)
  {
    allocate();
  }
  else
  {
    allocate(nrl,nrh,ncl,nch);
  }
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix::df1b2matrix(int nrl,int nrh)
{
  if (nrl>nrh)
  {
    allocate();
  }
  else
  {
    allocate(nrl,nrh);
  }
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix::df1b2matrix(int nrl,int nrh,const index_type& ncl,
  const index_type& nch)
{
  if (nrl>nrh)
  {
    allocate();
  }
  else
  {
    allocate(nrl,nrh,ncl,nch);
  }
}

/**
 * Description not yet available.
 * \param
 */
df1b2matrix::df1b2matrix(void)
{
  allocate();
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::allocate(int nrl,int nrh,int ncl,int nch, [[maybe_unused]] const char * s)
{
  allocate(nrl,nrh,ncl,nch);
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::allocate(int nrl,int nrh,const index_type& ncl,
  const index_type& nch, [[maybe_unused]] const char * s)
{
  allocate(nrl,nrh,ncl,nch);
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::allocate(int nrl,int nrh,int ncl,int nch)
{
  index_min=nrl;
  index_max=nrh;
  if ( (v = new df1b2vector[size()]) == 0)
  {
      cerr << " Error allocating memory in df1b2matrix contructor\n";
      ad_exit(21);
  }
  if ( (shape=new mat_shapex(v)) == 0)
  {
    cerr << " Error allocating memory in df1b2matrix contructor\n";
  }
  v -= indexmin();
  for (int i=nrl; i<=nrh; i++)
  {
    v[i].allocate(ncl,nch);
  }
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::allocate(int nrl,int nrh,const index_type& ncl,
  const index_type& nch)
{
  index_min=nrl;
  index_max=nrh;
  if ( (v = new df1b2vector[size()]) == 0)
  {
      cerr << " Error allocating memory in df1b2matrix contructor\n";
      ad_exit(21);
  }
  if ( (shape=new mat_shapex(v)) == 0)
  {
    cerr << " Error allocating memory in df1b2matrix contructor\n";
  }
  v -= indexmin();
  for (int i=nrl; i<=nrh; i++)
  {
    v[i].allocate(ncl(i),nch(i));
  }
}
/**
Copy constructor
*/
df1b2matrix::df1b2matrix(const df1b2matrix& x)
{
  index_min=x.index_min;
  index_max=x.index_max;
  v=x.v;
  shape=x.shape;
  if (shape) (shape->ncopies)++;
}
/**
Allocate square matrix.
\param nrl lower index
\param nrh upper index
*/
void df1b2matrix::allocate(int nrl,int nrh)
{
  index_min=nrl;
  index_max=nrh;
  if ( (v = new df1b2vector[size()]) == 0)
  {
      cerr << " Error allocating memory in df1b2matrix contructor\n";
      ad_exit(21);
  }
  if ( (shape=new mat_shapex(v)) == 0)
  {
    cerr << " Error allocating memory in df1b2matrix contructor\n";
  }
  v -= indexmin();
  /*
  for (int i=nrl; i<=nrh; i++)
  {
    v[i].allocate(ncl,nch);
  }
  */
}
/**
Destructor
*/
df1b2matrix::~df1b2matrix()
{
  if (shape)
  {
    if (shape->ncopies)
    {
      (shape->ncopies)--;
    }
    else
    {
      deallocate();
    }
  }
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::deallocate()
{
  if (shape)
  {
    v=(df1b2vector*)(shape->get_pointer());
    delete [] v;
    delete shape;
    allocate();
  }
}

/**
 * Description not yet available.
 * \param
 */
void df1b2matrix::allocate(void)
{
  index_min=1;
  index_max=0;
  v=0;
  shape=0;
}

/**
Copy vector values into column j of df1b2matrix.

\param j column index
\param values vector
*/
void df1b2matrix::colfill(const int j, const df1b2vector& values)
{
  //RETURN_ARRAYS_INCREMENT();   makes no sense for df1b2 stuff
  for (int i = rowmin(); i <= rowmax(); ++i)
  {
    (*this)[i][j] = values[i];
  }
  //RETURN_ARRAYS_DECREMENT();
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable sum(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=x(i);
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable mean(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  double fn=mmax-mmin+1;
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=x(i);
  }
  return tmp/fn;
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable norm2(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=square(x(i));
  }
  return tmp;
}
df1b2variable sumsq(const df1b2vector& _x) {return(norm2(_x));}

/**
 * Description not yet available.
 * \param
 */
df1b2variable norm(const df1b2vector& _x)
{
  ADUNCONST(df1b2vector,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=square(x(i));
  }
  return sqrt(tmp);
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable norm2(const df1b2matrix& _x)
{
  ADUNCONST(df1b2matrix,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=norm2(x(i));
  }
  return tmp;
}
df1b2variable sumsq(const df1b2matrix& _x) {return(norm2(_x));}

/**
 * Description not yet available.
 * \param
 */
df1b2variable norm(const df1b2matrix& _x)
{
  ADUNCONST(df1b2matrix,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=norm2(x(i));
  }
  return sqrt(tmp);
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable sum(const df1b2matrix& _x)
{
  ADUNCONST(df1b2matrix,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=sum(x(i));
  }
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
df1b2variable mean(const df1b2matrix& _x)
{
  ADUNCONST(df1b2matrix,x)
  df1b2variable tmp;
  tmp=0.0;
  int mmin=x.indexmin();
  int mmax=x.indexmax();
  int nitems=0;
  for (int i=mmin;i<=mmax;i++)
  {
    tmp+=sum(x(i));
    nitems+=x(i).indexmax()-x(i).indexmin()+1;
  }
#ifdef DEBUG
  assert(nitems > 0);
#endif
  return tmp/nitems;
}

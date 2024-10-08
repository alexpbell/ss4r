/**
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "fvar.h"
#include <string.h>
#include <ctype.h>
#ifdef DEBUG
  #include <cassert>
#endif
#include <climits>

/**
 * Description not yet available.
 * \param
 */
 struct dvec_ptr_ptr
 {
   void ** m;
 };

/**
Fill allocated dmatrix with values from input parameter s.

    dmatrix m(1, 2, 1, 3);
    char values[] = "{0, 1, 2} {3, 4, 5}";
    m.fill(values);

\param s values
*/
void dmatrix::fill(const char* s)
{
  #ifdef DIAG
    myheapcheck("Entering dmatrix(const char * s)" );
  #endif

  size_t len = strlen(s);
#ifdef DEBUG
  assert(len <= INT_MAX);
#endif
  int n = static_cast<int>(len);
  int braces = 0;
  int nrow = 0;
  int ncol = 0;
  int mincolumn=colmin();
  int minrow=rowmin();
  rowshift(1);
  colshift(1);

  const int MAXROWS = 5050;
  ivector columns(1, MAXROWS);
  ivector k1(1, MAXROWS);
  ivector k2(1, MAXROWS);

  for (int k = 0; k < n; k++)
  {
    if (s[k] == '{')
    {
      braces ++;
      if (braces != 1)
      {
        cerr << "Unbalanced braces in dmatrix::fill(const char * s)\n";
        cerr << "at character " << k << "\n";
        ad_exit(1);
      }
      ncol = 1;
      k1[nrow+1] = k;
    }
    else if (s[k] == '}')
    {
      braces --;
      if (braces != 0)
      {
        cerr << "Unbalanced braces in dmatrix::dmatrix(const char * s)\n";
        cerr << "at character " << k << "\n";
        ad_exit(1);
      }
      k2[nrow+1] = k;
      nrow ++;
      if (nrow > MAXROWS)
      {
        cerr << "Too many rows in dmatrix::dmatrix(const char * s)\n";
        ad_exit(1);
      }
      columns[nrow] = ncol;
    }
    else if (s[k] == ',')
    {
      if (braces != 0)
      {
        ncol++;
      }
    }
  }

  if (braces != 0)
  {
    cerr << "Unbalanced braces in dmatrix::dmatrix(const char * s)\n";
    cerr << s << "\n";
    ad_exit(1);
  }

  {
    unsigned int _nrow = static_cast<unsigned int>(nrow);
    unsigned int _rowsize = rowsize();
    if (_nrow < _rowsize)
    {
      cerr << " Not enough rows in the data for dmatrix::fill(const char *) \n";
      ad_exit(1);
    }
    else if (_nrow > _rowsize)
    {
      cerr << " Too many rows in the data for dmatrix::fill(const char *) \n";
      ad_exit(1);
    }
  }

  for (int i=1; i<=nrow; i++)
  {
    int index = rowmin() + i - 1;
    unsigned int size = ((*this)[index]).size();
    unsigned int cols = static_cast<unsigned int>(columns[i]);
    cout << "row " << i << " matrix " << size << " colvector " << cols << "\n";

    if (size > cols)
    {
      cerr << " Not enough columns in the data in row "
           << i << " for dmatrix::fill(const char *) \n";
      ad_exit(1);
    }
    else if (size < cols)
    {
      cerr << " Too many columns in the data in row "
           << i << " for dmatrix::fill(const char *) \n";
      ad_exit(1);
    }
  }

  for (int i=rowmin(); i<=rowmax(); i++)
  {
    char *t = new char[strlen(s)+1];
    //t = (char*) new[strlen(s)+1];

    for (int k = k1[i]; k <= k2[i]; k++)
    {
      t[k-k1[i]] = s[k];
    }
    t[k2[i]-k1[i]+1] = '\0';

    dvector tt(t);
    tt.shift(( (*this)(i)).indexmin() );

   //    (*this)(i)=dvector(tt);
    (*this)(i)=tt;

    delete[] t;
    t = 0;
  }
  rowshift(minrow);
  colshift(mincolumn);
}

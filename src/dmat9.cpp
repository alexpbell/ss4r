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
#endif

/**
 * Description not yet available.
 * \param
 */
mat_shape::mat_shape(int rl,int ru,int cl,int cu)
{
#ifdef DEBUG
  assert(ru >= rl);
  assert(cu >= cl);
#endif
  row_min=rl;
  row_max=ru;
  col_min=cl;
  col_max=cu;
  nrows=(unsigned int)(ru-rl+1);
  ncols=(unsigned int)(cu-cl+1);
  ncopies=0;
}

/**
Changes the range of valid indices for the rows.
\param min shifts column values index_min and index_max.
*/
void mat_shape::colshift(int min)
{
  col_max = col_max - col_min + min;
  col_min = min;
}
/**
Changes the range of valid indices for the rows.
\param min shifts row values index_min and index_max.
*/
void dmatrix::rowshift(int min)
{
  m = m + rowmin() - min;
  index_max += min - index_min;
  index_min = min;
}
/**
 * Description not yet available.
 * \param
 */
void mat_shape::rowshift(int min)
{
  row_max=row_max-row_min+min;
  row_min=min;
}

/**
 * Description not yet available.
 * \param
 */
 void dmatrix::colshift( int min)
 {
   for (int i=rowmin(); i<=rowmax(); i++)
   {
     this->elem(i).shift(min);
   }
   //shape->colshift(min);
 }

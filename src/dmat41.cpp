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
void banded_symmetric_dmatrix::initialize(void)
{
  for (int i=rowmin();i<=rowmax();i++)
  {
    (*this)(i).initialize();
  }
}

/**
 * Description not yet available.
 * \param
 */
void banded_lower_triangular_dmatrix::initialize(void)
{
  for (int i=rowmin();i<=rowmax();i++)
  {
    (*this)(i).initialize();
  }
}

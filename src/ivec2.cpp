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
void ivector::initialize(void)
{
  int min = indexmin();
  int max = indexmax();
  int* pvi = v + min;
  for (int i = min; i <= max; ++i)
  {
    *pvi = 0;
    ++pvi;
  }
}

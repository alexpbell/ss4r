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
//#undef OPT_LIB
#include "fvar.h"

#if !defined(OPT_LIB)

/**
 * Description not yet available.
 * \param
 */
 dvar_vector& dvar_matrix::operator() (int i)
 {
     if (i<rowmin())
     {
       cerr << "matrix bound exceeded -- row index too low in "
       "dvar_matrix::operator()" << "value was" << i << endl;
       ad_exit(21);
     }
     if (i>rowmax())
     {
       cerr << "matrix bound exceeded -- row index too high in "
       "dvar_matrix::operator()" << "value was" << i << endl;
       ad_exit(22);
     }
   return m[i];
 }
#endif

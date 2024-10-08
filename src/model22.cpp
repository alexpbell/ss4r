/**
@file
@author David Fournier
@copyright Copyright (c) 2008-2020 Regents of the University of California

@brief Function named_dvar3_array::allocate
*/
#include "admodel.h"

void named_dvar3_array::allocate(const ad_integer& sl,
  const ad_integer& sh,const index_type& nrl,const index_type& nrh,
  const index_type& ncl,const index_type& nch,const char * s)
{
  dvar3_array::allocate(sl,sh,nrl,nrh,ncl,nch);
  model_name_tag::allocate(s);
}

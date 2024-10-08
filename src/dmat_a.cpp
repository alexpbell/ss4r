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
 struct dvec_ptr_ptr
 {
   void ** m;
 };

/**
 * Description not yet available.
 * \param
 */
 dmatrix::dmatrix(const dvar_matrix_position& pos)
 {
   int nrl=pos.row_min;
   int nrh=pos.row_max;
   const ivector& ncl=pos.lb;
   const ivector& nch=pos.ub;
#ifdef DIAG
   myheapcheck("Entering dmatrix(nrl,nrh,ncl,nch)" );
#endif

#ifndef OPT_LIB
   if (nrl != ncl.indexmin() || nrh != ncl.indexmax()
     || nrl != nch.indexmin() || nrh != nch.indexmax())
   {
     cerr << "Incompatible array bounds in "
     "dmatrix(int nrl,int nrh, const ivector& ncl, const ivector& nch)\n";
     ad_exit(1);
   }
#endif
   index_min=nrl;
   index_max=nrh;

   if ( (m = new dvector[rowsize()]) == 0)
   {
     cerr << " Error allocating memory in dmatrix contructor\n";
     ad_exit(21);
   }
   if ( (shape = new mat_shapex(m)) == 0)
   {
     cerr << " Error allocating memory in dmatrix contructor\n";
     ad_exit(21);
   }

#ifdef DIAG
   cerr << "Created a dmatrix with adress "<< farptr_tolong(m)<<"\n";
#endif

   m -= rowmin();

   dvector* pm = m + nrl;
   int* pncl = ncl.get_v() + nrl;
   int* pnch = nch.get_v() + nrl;
   for (int i = nrl; i <= nrh; ++i)
   {
     pm->allocate(*pncl, *pnch);
     ++pm;
     ++pncl;
     ++pnch;
#ifdef DIAG
     cerr << "Created a dvector with address "<< farptr_tolong((void*)(m+i))<<"\n";
#endif
   }

#ifdef DIAG
   myheapcheck("Leaving dmatrix(nrl,nrh,ncl,nch)" );
#endif
 }

/**
 * Description not yet available.
 * \param
 */
 dmatrix::dmatrix(const dmatrix_position& pos)
 {
   int nrl=pos.row_min;
   int nrh=pos.row_max;
   const ivector& ncl=pos.lb;
   const ivector& nch=pos.ub;

#ifdef DIAG
     myheapcheck("Entering dmatrix(nrl,nrh,ncl,nch)" );
#endif

#ifndef OPT_LIB
   if (nrl != ncl.indexmin() || nrh != ncl.indexmax()
       || nrl != nch.indexmin() || nrh != nch.indexmax())
   {
     cerr << "Incompatible array bounds in "
     "dmatrix(int nrl,int nrh, const ivector& ncl, const ivector& nch)\n";
     ad_exit(1);
   }
#endif
   index_min=nrl;
   index_max=nrh;

   if ( (m = new dvector[rowsize()]) == 0)
   {
     cerr << " Error allocating memory in dmatrix contructor\n";
     ad_exit(21);
   }

   if ( (shape = new mat_shapex(m))== 0)
   {
     cerr << " Error allocating memory in dmatrix contructor\n";
     ad_exit(21);
   }

#ifdef DIAG
     cerr << "Created a dmatrix with adress "<< farptr_tolong(m)<<"\n";
#endif

   m -= rowmin();

   dvector* pm = m + nrl;
   int* pncl = ncl.get_v() + nrl;
   int* pnch = nch.get_v() + nrl;
   for (int i = nrl; i <= nrh; ++i)
   {
     pm->allocate(*pncl, *pnch);
     ++pncl;
     ++pnch;
     ++pm;
   }
#ifdef DIAG
     myheapcheck("Leaving dmatrix(nrl,nrh,ncl,nch)" );
#endif
 }

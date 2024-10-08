/**
\file

Author: David Fournier
Copyright (c) 2008-2012 Regents of the University of California
*/
#include "fvar.h"

#include <stdlib.h>
#ifdef DEBUG
  #include <cassert>
#endif

/**
Default constructor
*/
dlist::dlist()
{
  int on,nopt = 0;
  if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-mdl",nopt))>-1)
  {
    if (nopt == 1)
    {
      int i = atoi(ad_comm::argv[on+1]);
      if (i > 0)
        gradient_structure::MAX_DLINKS = (unsigned int)i;
    }
    else
    {
      cerr << "Wrong number of options to -mdl -- must be 1"
        " you have " << nopt << endl;
      ad_exit(1);
    }
  }

  last = NULL;
  nlinks = 0;
  dlink_addresses = new dlink*[gradient_structure::MAX_DLINKS];
  ddlist_space =
    (char*)malloc(gradient_structure::MAX_DLINKS * sizeof(dlink));
  variables_save = new double[gradient_structure::MAX_DLINKS];

#ifdef DEBUG
  //fails for insufficient memory to allocate space for dvariables save buffer
  assert(variables_save != NULL);
#endif

  //Initialize addresses to zero
  memset(dlink_addresses, 0, sizeof(dlink*) * gradient_structure::MAX_DLINKS);
}
/**
Destructor
*/
dlist::~dlist()
{
  if (dlink_addresses)
  {
    delete [] dlink_addresses;
    dlink_addresses = NULL;
  }
  if (ddlist_space)
  {
    ::free(ddlist_space);
    ddlist_space = NULL;
  }
  if (variables_save)
  {
    delete [] variables_save;
    variables_save = NULL;
  }
}
/**
Return new unlinked node.
*/
dlink* dlist::create()
{
#ifdef DEBUG
  //If fails, then need to increase the maximum number of dlinks.
  assert(nlinks < gradient_structure::MAX_DLINKS);
#endif

  dlink* link = (dlink*)(&ddlist_space[sizeof(dlink) * nlinks]);

#ifdef DEBUG
  assert(link);
#endif

  //Do not add to list.
  link->prev = NULL;

  //Keep track of the links so you can zero them out (ie gradcalc).
  dlink_addresses[nlinks] = link;
  ++nlinks;

  return link;
}
/**
If list is not empty, pop and return last node.

\return 0 empty list.
*/
dlink* dlist::last_remove()
{
  dlink* link = last;
  if (link)
  {
    last = link->prev;
    link->prev = NULL;
  }
  return link;
}
/**
Append link to list.

\param link node
*/
dlink* dlist::append(dlink* link)
{
#ifdef DEBUG
  //Should fail if link is NULL.
  assert(link);
#endif

  link->prev = last;
  last = link;

  return last;
}
void dlist::initialize()
{
  dlink** dest = dlink_addresses;
  for (unsigned int i = 0; i < nlinks; ++i)
  {
    (*dest)->di.x = 0;
    ++dest;
  }
}
/**
Save variables to a buffer.
*/
void dlist::save_variables()
{
  dlink** src = dlink_addresses;
  double* dest = variables_save;
  for (unsigned int i = 0; i < nlinks; ++i)
  {
    *dest = (*src)->di.x;
    ++dest;
    ++src;
  }
}
/**
Restore variables from buffer.
*/
void dlist::restore_variables()
{
  dlink** dest = dlink_addresses;
  double* src = variables_save;
  for (unsigned int i = 0; i < nlinks; ++i)
  {
    (*dest)->di.x = *src;
    ++dest;
    ++src;
  }
}
/**
Get total addresses stored.
*/
size_t dlist::total_addresses() const
{
  size_t total = 0;
  for (unsigned int i = 0; i < gradient_structure::MAX_DLINKS; ++i)
  {
    if (dlink_addresses[i] != 0)
    {
      total++;
    }
  }
  return total;
}
/**
Check link list integrity.
*/
void dlist::check_list(void)
{
  dlink* tmp_last=last;

  unsigned int count=0;
  while(tmp_last && count <=nlinks)
  {
    count+=1;
    if (count > nlinks)
    {
      cerr << "In check_list() number of links > number created\n";
      cerr << " The number created was "<< nlinks << endl;
    }

    dlink* tmp = tmp_last->prev;

//  cout << "last =" << _farptr_tolong(last) << "\n";
//  cout << "last->prev =" << _farptr_tolong(last->prev) << "\n";
//  cout << "deleted dlink with address" << _farptr_tolong(last) << "\n";

    tmp_last = tmp;
  }
  cerr << "In check_list() number of free links is " << count << endl;
}

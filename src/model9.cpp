/*
 * $Id$
 *
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "fvar.h"
#include "admodel.h"

int operator+(const int n, const data_int& v)
{
  return n+v.val;
}

int operator+(const data_int& v, const int n)
{
  return n+v.val;
}

int operator+(const data_int& v, const data_int& n)
{
  return n.val + v.val;
}

data_int& data_int::operator=(const int xx)
{
  val=xx;
  return *this;
}

ad_integer::ad_integer(const data_int& _d) : d(int(*(data_int*)(&_d))) {}
/*

index_type::index_type(const data_int& _x)
{
  p = new number_index(int((data_int&)(_x)));
}
*/

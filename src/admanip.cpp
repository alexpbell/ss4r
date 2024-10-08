/**
@file
@author David Fournier
@copyright Copyright (c) 2008-2020 Regents of the University of California

@brief Overloaded ostream outputs functions.
*/
#include "fvar.h"

/*
ostream& setfixed(const ostream& _s)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::fixed,ios::floatfield);
  return s;
}

ostream& setscientific(const ostream& _s)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::scientific,ios::floatfield);
  return s;
}


ostream& setshowpoint(const ostream& _s)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::showpoint);
  return s;
}
*/

/**
 * Description not yet available.
 * \param
 */
preshowpoint setshowpoint(void)
{
  preshowpoint tmp;
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
ostream& operator<<(const ostream& _s, [[maybe_unused]] preshowpoint p)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::showpoint);
  return s;
}

/**
 * Description not yet available.
 * \param
 */
prefixed setfixed(void)
{
  prefixed tmp;
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
ostream& operator<<(const ostream& _s, [[maybe_unused]] prefixed p)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::fixed,ios::floatfield);
  return s;
}

/**
 * Description not yet available.
 * \param
 */
prescientific setscientific(void)
{
  prescientific tmp;
  return tmp;
}

/**
 * Description not yet available.
 * \param
 */
ostream& operator<<(const ostream& _s, [[maybe_unused]] prescientific p)
{
  ostream& s=(ostream&)(_s);
  s.setf(ios::scientific,ios::floatfield);
  return s;
}

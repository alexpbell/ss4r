/**
@file
@author David Fournier
@copyright Copyright (c) 2008-2020 Regents of the University of California

@brief Functions for prevariable to compute cube(const prevariable&) and fourth(const prevariable&).
*/

#include "fvar.h"

/**
Returns variable result of v1 cubed.

\ingroup misc
\param v1 variable base
\return \f$v^3\f$
*/
prevariable& cube(const prevariable& v1)
{
  if (++gradient_structure::_instance->RETURN_PTR > gradient_structure::_instance->MAX_RETURN) 
    gradient_structure::_instance->RETURN_PTR = gradient_structure::_instance->MIN_RETURN;

  double x=value(v1);
  double x2=x*x;

  gradient_structure::_instance->RETURN_PTR->v->x=x2*x;
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation2,
    &(gradient_structure::_instance->RETURN_PTR->v->x), &(v1.v->x), 3.0*x2);

  return *gradient_structure::_instance->RETURN_PTR;
}
/**
Returns variable result of v1 raised to the power of four.

\ingroup misc
\param v1 variable base
\return \f$v^4\f$
*/
prevariable& fourth(const prevariable& v1)
{
  if (++gradient_structure::_instance->RETURN_PTR > gradient_structure::_instance->MAX_RETURN) 
    gradient_structure::_instance->RETURN_PTR = gradient_structure::_instance->MIN_RETURN;

  double x=value(v1);
  double x2=x*x;

  gradient_structure::_instance->RETURN_PTR->v->x=x2*x2;
  gradient_structure::GRAD_STACK1->set_gradient_stack(default_evaluation2,
    &(gradient_structure::_instance->RETURN_PTR->v->x), &(v1.v->x), 4.0*x2*x);

  return *gradient_structure::_instance->RETURN_PTR;
}
